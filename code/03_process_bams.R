library(rbamtools)
library(parallel)
library(lme4)
library(stringdist)
library(mgcv)
library(bio3d)
library(msm)
library(foreach)
library(doMC)
registerDoMC(cores=30)

t.var=c(0,24,48,72,96)
out.base.dir='/media/jbloom/d1/coupled_CRISPR/Experiments/082316/'

readID=list.files(paste0(out.base.dir, 'fastq/'))
# subsetting specific to this experiment
readID=unique(gsub('_R.*.fq.gz', '', readID))[1:12]

# output here 
dir.create(paste(out.base.dir, 'processed/RData/', sep=''))

source('/media/jbloom/d1/coupled_CRISPR/code/accessory_functions.R')

## Extract all relevant information from BAM files and restructure as R Data files ---------------------------------
# note that we structure as a very big list here, but this can be a problem if output from hiseq is huge
## also build a crude list of counts per oligo (will refine later) 
all.counts=list()
for( rr in readID) {
    print(rr)
    bam.file=paste(out.base.dir, 'processed/bam/', rr, '.bam', sep='')

    # use rbamtools to extract aligned results 
    reader = bamReader(bam.file,idx=TRUE)
    contigs = (refSeqDict(getHeaderText(reader)))@SN
    contigs.length=(refSeqDict(getHeaderText(reader)))@LN

    # construct a table to parse the information about each oligo construct
    ac=bamCountAll(reader)
    bamClose(reader)

    reader=bamReader(bam.file,idx=TRUE)
    branges=list()
    pb =txtProgressBar(min = 1, max =length(contigs), style = 3)
    for(n in 1:length(contigs)) {
        setTxtProgressBar(pb, n)
        # argh, 0 indexing
        branges[[contigs[n]]]=bamRange(reader, coords=c(n-1, 1, contigs.length[n])) 
    }
    close(pb)
    bamClose(reader)

    bdf=mclapply(branges, function(x) {
           a=as.data.frame(x)   
           # modified 062916
           # modified 082416 for hiseq
           a$barcode=sapply(strsplit(a$name, ':'), function(x)x[2])
           a$barcode.type=sapply(strsplit(a$name, ':'), function(x)x[3])
           return(a)
    }, mc.cores=66)
    all.counts[[rr]]=ac
    saveRDS(bdf, file=paste0(out.base.dir, 'processed/RData/', rr, '.RDS'))
}   
#------------------------------------------------------------------------------------------------------------------
 
#save(all.counts, file=paste0(out.base.dir, '/processed/RData/all.counts.RData') )
load( paste0(out.base.dir, '/processed/RData/all.counts.RData') )

#report number of 'expected' oligos with non-zero counts
sapply(all.counts, function(x) sum(x$nAligns>0))

#total counts per sample
sapply(all.counts, function(x) sum(x$nAligns))

# load all the raw data and create one very large R list structure  ----------------------------------
expt=list()
for(rr in readID) {  print(rr)
                     expt[[rr]]=readRDS(paste0(out.base.dir, 'processed/RData/', rr, '.RDS'))    }
names(expt)=readID
names(expt)=gsub('72-plate', '72plate', names(expt))
#saveRDS(expt, file = paste0(out.base.dir, 'processed/RData/expt.RDS'))
#-----------------------------------------------------------------------------------------------------

expt=readRDS(file=paste0(out.base.dir, 'processed/RData/expt.RDS'))


# Each lane (element of expt) is a mixture of two experiments. This block of code will split reads based on the internal barcode sequence into separate data structures
    # Key for experiment names M. Sadhu 082316
    #106 is WT 
    #131 is NMD-
    #295 = swws
    #296 = wssw
   # swws=mclapply(expt, function(x){ lapply(x, function(y) { return(y[y$barcode.type=='swws',])})}, mc.cores=66)
    swws=lapply(expt, function(x){ lapply(x, function(y) { return(y[y$barcode.type=='swws',])})})
    
    #wssw=mclapply(expt, function(x){ lapply(x, function(y) { return(y[y$barcode.type=='wssw',])})}, mc.cores=66)
    wssw=lapply(expt, function(x){ lapply(x, function(y) { return(y[y$barcode.type=='wssw',])})})

    expt.key=do.call('rbind', strsplit(names(expt), ':|_|-'))
    expt.key[expt.key=='106'] = 'WT'
    expt.key[expt.key=='131'] = 'NMDm'
    expt.key[expt.key=='295'] = 'swws'
    expt.key[expt.key=='296'] = 'wssw'

    # rename experiment to WT or NMD and time point
    swws.key=apply(expt.key, 1, function(x) paste( x[1],'_', x[6], sep=''))
    wssw.key=apply(expt.key, 1, function(x) paste( x[3],'_', x[6], sep=''))

    names(swws)=swws.key
    names(wssw)=wssw.key

    swws=swws[sort(names(swws))]
    wssw=wssw[sort(names(wssw))]

# Now combine experiments back together ------------------------------------------------------------------------------------------------------------------------------

c.expt=list()
for(i in 1:length(swws)) {
    print(i)
    c.expt[[names(swws)[i]]]=mapply('rbind', swws[[i]], wssw[[i]], SIMPLIFY=FALSE)
}
#saveRDS(c.expt, file = paste0(out.base.dir, 'processed/RData/combined_expt.RDS'))

# read in list structure
c.expt=readRDS(file = paste0(out.base.dir, 'processed/RData/combined_expt.RDS'))

# this code will collapse barcodes per oligo by edit distance (if within 3 bases then barcode is likely the same ... look at this distribution)
# also error correct sequence by barcode, by some criteria if the predominate species for a given barcode is one sequence then deviations from that sequence are likely errors
#i=300
#54 and 22 
#'ATAGGGTCGTAC 10
#'TGTCAGTCGTAC 33'
#check out oligo 105 barcode #32
barcode.grouping=mclapply(1:nrow(oligos),function(i){
    #for(i in 1:nrow(oligos)) {
    print(i)
    u1=unlist(sapply(c.expt, function(x) x[[i]]$barcode))
    if(is.null(u1) | length(u1)==0){
        print('zero barcodes')
       # barcode.grouping[[i]]=
            return(data.frame(barcode='', barcode.count=0, bgroup=0, stringsAsFactors=F))
        #next;
    } else {
        ubc=rle(sort(as.vector(u1)))
    }
    if(length(ubc$lengths)==1) {
      #  barcode.grouping[[i]]=
            return(data.frame(barcode=ubc$values[1], barcode.count=ubc$lengths[1], bgroup=1, stringsAsFactors=F))
    } else {
        hc= hclust(as.dist(adist(ubc$values)))
       
        b.groups=cutree(hc, h=3)
        #rbgs=rle(sort(b.groups))
        #rbgsgt1=rbgs$values[which(rbgs$lengths>1)]
        #rect.hclust(hc, h=3, which=rbgsgt1)
        ad=data.frame(barcode=ubc[[2]], barcode.count=ubc[[1]], bgroup=b.groups, stringsAsFactors=F)
        
        #nname=paste0('[', sprintf("%3d", ad$bgroup), '] ', ubc$values, sprintf("%4d",ad$barcode.count))
        #par(family='mono')
        #plot(hc, labels=nname, xlab='hierarchical clustering of barcodes based on edit distance', cex=1.2, main=paste('oligo #32 of 10971'))

        barcode.grouping[[i]]=
            return(ad[order(ad$bgroup, ad$barcode.count),])
    }
}, mc.cores=10)
names(barcode.grouping)=names(c.expt[[1]])


#  !! call by reference only 
for(i in 1:nrow(oligos)) {
    print(i)
    bg=barcode.grouping[[i]]
    for(j in 1:length(c.expt)) {
        bg.vec= bg$bgroup[match(c.expt[[j]][[i]]$barcode, bg$barcode)]
        if(length(bg.vec)==0) {
             c.expt[[j]][[i]]$barcode.group=list()
        }    else {
             c.expt[[j]][[i]]$barcode.group=bg.vec 
        }
    }
}   
#saveRDS(c.expt, file = paste0(out.base.dir, 'processed/RData/combined_expt_bc_error_corrected.RDS'))

c.expt=readRDS(file = paste0(out.base.dir, 'processed/RData/combined_expt_bc_error_corrected.RDS'))

# record the  counts of unique sequences for each of the unique barcode sequences 
seq.tables=list()
#seq.tables=mclapply( 1:nrow(oligos), function(i) {
for(i in 1:nrow(oligos)) {
    print(i)
    #now error correct reads
    all.seqs=unlist(sapply(c.expt, function(x) x[[i]]$seq))
    all.qual=unlist(sapply(c.expt, function(x) x[[i]]$qual))
    all.cigar=unlist(sapply(c.expt, function(x) x[[i]]$cigar))
    all.revstrand=unlist(sapply(c.expt, function(x) x[[i]]$revstrand))
    all.barcode=unlist(sapply(c.expt, function(x) x[[i]]$barcode))
    all.barcodetype=unlist(sapply(c.expt, function(x) x[[i]]$barcode.type))

    all.bcgroups=unlist(sapply(c.expt, function(x) x[[i]]$barcode.group))
    all.experiment=rep(names(c.expt), sapply(c.expt, function(x) nrow(x[[i]])))
    asb=data.frame(seqs=as.vector(all.seqs), 
                   qual=as.vector(all.qual),
                   cigar=as.vector(all.cigar),
                   revstrands=as.vector(all.revstrand),
                   barcode.group=as.vector(all.bcgroups),
                   barcode.type=as.vector(all.barcodetype),
                   barcode=as.vector(all.barcode),
                   experiment=as.vector(all.experiment)   ,           
               stringsAsFactors=F)
    if(nrow(asb)==0) {
        #return(
        seq.tables[[names(c.expt[[1]])[i]]]=
                   data.frame(
                                                seqs         =NULL,  
                                                qual         =NULL,   
                                                cigar        =NULL,   
                                                revstrands   =NULL,   
                                                barcode.group  =NULL, 
                                                barcode.type   =NULL, 
                                                barcode        =NULL, 
                                                max.obs.seq.cnt=NULL, 
                                                WT_0           =NULL, 
                                                WT_24          =NULL, 
                                                WT_48          =NULL, 
                                                WT_72          =NULL, 
                                                WT_96          =NULL, 
                                                WT_72plate     =NULL, 
                                                NMDm_0         =NULL, 
                                                NMDm_24        =NULL, 
                                                NMDm_48        =NULL, 
                                                NMDm_72        =NULL, 
                                                NMDm_96        =NULL, 
                                                NMDm_72plate   =NULL) 
                   #)
    } else {

        asb.split=split(asb, list(asb$barcode.group))
        asb.max.seq=(lapply(asb.split, 
                   function(y) { sort(table(y$seqs), decreasing=T)[1] } ))
        ams=sapply(asb.max.seq, function(x) x[1])
        names(ams)=as.vector(sapply(asb.max.seq, names))
        
        best.seq=do.call('rbind', lapply(1:length(ams), function(k) {
            kout=asb.split[[k]][match(names(ams)[k], asb.split[[k]]$seqs),]
            kout$max.obs.seq.cnt=as.vector(ams[k])
            return(kout)
        } ))
        rownames(best.seq)=NULL
        best.seq$experiment=NULL
        
        tableme=table(asb$barcode.group, asb$experiment)
        table.out=matrix(0, nrow(best.seq), length(names(c.expt)))
        colnames(table.out)=names(c.expt)[c(7:10,12,11,1:4,6,5)]
        cmi=match(colnames(tableme), colnames(table.out))
        table.out[,cmi]=tableme
        #return(
        seq.tables[[names(c.expt[[1]])[i]]]=
               data.frame(best.seq, table.out, stringsAsFactors=F)
               #print(seq.tables[[names(c.expt[[1]])[i]]])
           #)
    }
}
#saveRDS(seq.tables,  file = paste0(out.base.dir, 'processed/RData/seq.tables.RData'))


giant.table=do.call('rbind', seq.tables)
giant.table$oligo=rep(names(seq.tables), sapply(seq.tables, nrow))
#saveRDS(giant.table,  file = paste0(out.base.dir, 'processed/RData/giant.table.RData'))
