library(rbamtools)
library(parallel)
library(lme4)
library(stringdist)
library(mgcv)
library(bio3d)
library(msm)
library(Biostrings)
library(seqinr)
library(yeastExpData)
library(gdata)
library("BSgenome.Scerevisiae.UCSC.sacCer3")
sacCer3=BSgenome.Scerevisiae.UCSC.sacCer3
library("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene")
txdb=TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
library(VariantAnnotation)
library(foreach)
library(doMC)
registerDoMC(cores=30)
library(data.table)

t.var=c(0,24,48,72,96)
out.base.dir='/media/jbloom/d1/coupled_CRISPR/Experiments/082316/'

readID=list.files(paste0(out.base.dir, 'fastq/'))
# subsetting specific to this experiment
readID=unique(gsub('_R.*.fq.gz', '', readID))[1:12]

# output here 
dir.create(paste(out.base.dir, 'processed/RData/', sep=''))

#split out functions necessary from pre-processing (03_ vs 04_) 
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

c.expt =list()
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
        #investigate using stringdist here Frank Albert 12/15/16
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
        
        # this grabs all info for most common barcode 
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

seq.tables=readRDS(file = paste0(out.base.dir, 'processed/RData/seq.tables.RData'))
giant.table=readRDS(file = paste0(out.base.dir, 'processed/RData/giant.table.RData'))
# this is convenient
rownames(oligos)=names(seq.tables)

giant.table$start.align=0                    # start.align=start(subject(pas)),
giant.table$end.align=0                      # end.align=end(subject(pas)),
giant.table$n.indel=0                        # n.indel=indel.cnt$counts,
giant.table$n.insertion=0                    # n.insertion=insertion.cnt$counts,
giant.table$n.mismatch=0                     # n.mismatch=nmismatch(pas),
giant.table$indel.pos=''                     # indel.pos=indel.pos,
giant.table$insertion.pos=''                 # insertion.pos=insertion.pos,
giant.table$mm.pos=''                        # mm.pos=mm.pos, 
giant.table$pid=0                            # pid=pid(pas),
giant.table$pam.indel.cnt=0                  # pam.indel.cnt=pam.indel.cnt,
giant.table$pam.mm.cnt=0                     # pam.mm.cnt= pam.mm.cnt, 
giant.table$u1.indel.cnt=                    # u1.indel.cnt=u1.indel.cnt, 
giant.table$u1.mm.cnt=0                      # u1.mm.cnt=u1.mm.cnt, 
giant.table$u.indel.cnt=0                    # u.indel.cnt=u.indel.cnt, 
giant.table$u.mm.cnt=0                       # u.mm.cnt=u.mm.cnt, 
giant.table$d.indel.cnt=0                    # d.indel.cnt=d.indel.cnt, 
giant.table$d.mm.cnt=0                       # d.mm.cnt=d.mm.cnt, 
giant.table$d1.indel.cnt=0                   # d1.indel.cnt=d1.indel.cnt, 
giant.table$d1.mm.cnt=0                      # d1.mm.cnt=d1.mm.cnt,
giant.table$alignment=''                     # alignment=as.character(apas))


#giant.table.perfect=giant.table$cigar=='101M'
#giant.table.whack =giant.table$cigar!='101M'
#gt=giant.table[!giant.table.whack,]
#gtw=giant.table[giant.table.whack,]
gtwo=split(giant.table, giant.table$oligo)

#o1=names(gtwo)[2]
window.size=20
alength=101
pam.window.ind=c(50:52)
u.window.ind=c(c((min(pam.window.ind)-1)-window.size) : c(min(pam.window.ind)-1))
d.window.ind=c( (max(pam.window.ind)+1) : ((max(pam.window.ind)+1)+window.size) )
u1.window.ind=c(1 : ( min(u.window.ind)-1))
d1.window.ind=c(c(max(d.window.ind)+1) : alength)

registerDoMC(cores=70)
#realigned.reads=list()
#for(o1 in names(gtwo)) {

#o1=names(which.max(sapply(gtwo, nrow)))
realigned.reads=foreach(o1=names(gtwo)) %dopar% {
#realigned.reads=list()
#for(o1 in names(gtwo)) {
    print(o1)
    oset=gtwo[[o1]]
   
    subject.seq=oligos[o1,]$repair.coding.strand
    subject.coding.strand=subject.seq==oligos[o1,]$stopPAM_oligos
        
    # patternQuality=SolexaQuality(oset$qual)
    observed.seqs=oset$seqs
    fobserved.seqs=as.factor(observed.seqs)
    to.align.seqs=as.character(levels(fobserved.seqs))
    
    mat = nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = FALSE, type='DNA')
          
    if(subject.coding.strand) {
                  pas=pairwiseAlignment(pattern=DNAStringSet(to.align.seqs), 
                                        subject=subject.seq,type='local', substitutionMatrix=mat, gapOpening=5, gapExtension=2)
    } else {
        pas=pairwiseAlignment(pattern=reverseComplement(DNAStringSet(to.align.seqs)), 
                              subject=subject.seq, type='local',substitutionMatrix=mat, gapOpening=5, gapExtension=2)
    }
    apas=aligned(pas)
  
    apas.mat=as.matrix(apas)
    mismatches=mismatchTable(pas)
    #msplit=split(mismatches$PatternStart, mismatches$PatternId)
    msplit=split(mismatches$SubjectStart, mismatches$PatternId)
    
    mm.cnt=list()
    mm.cnt$pos=vector("list", length(pas))
    names(mm.cnt$pos)=as.character(seq(1:length(pas)))
    mm.cnt$pos[names(msplit)]=msplit
    mm.cnt$counts=as.vector(sapply(mm.cnt$pos, length))
    
    deletion.intervals=RangedData(deletion(pas))

    subject.start=start(subject(pas))
    deletion.intervals=RangedData(deletion(pas))
    deletion.intervals$subject.start=subject.start[as.numeric(deletion.intervals$space)]
    deletion.intervals$start=start(deletion.intervals)+deletion.intervals$subject.start-1
    deletion.intervals$end= end(deletion.intervals)+deletion.intervals$subject.start-1
    sdi=as.character(space(deletion.intervals)) 
    sds=as.numeric(deletion.intervals$start)
    eds=as.numeric(deletion.intervals$end)

    insertion.intervals=RangedData(insertion(pas))
    insertion.intervals$subject.start=subject.start[as.numeric(insertion.intervals$space)]
    insertion.intervals$start=start(insertion.intervals)+insertion.intervals$subject.start-1
    insertion.intervals$end= end(insertion.intervals)+insertion.intervals$subject.start-1
    sii=as.character(space(insertion.intervals)) 
    sis=as.numeric(insertion.intervals$start)
    eis=as.numeric(insertion.intervals$end)


    #insertion.intervals=RangedData(insertion(pas))
    
    indel.cnt=list()
    indel.cnt$pos=vector("list", length(pas))
    names(indel.cnt$pos)=as.character(seq(1:length(pas)))

    insertion.cnt=list()
    insertion.cnt$pos=vector("list", length(pas))
    names(insertion.cnt$pos)=as.character(seq(1:length(pas)))

    if(nrow(deletion.intervals)>0) {
        for(i in 1:length(sdi)) {
            indel.cnt$pos[[sdi[i]]]=c(indel.cnt$pos[[sdi[i]]], c(sds[i]:eds[i]))
        }
    }

    if(nrow(insertion.intervals)>0) {
        for(i in 1:length(sii)) {
            indel.cnt$pos[[sii[i]]]=c(indel.cnt$pos[[sii[i]]], c(sis[i]:eis[i]))
        }
        for(i in 1:length(sii)) {
            insertion.cnt$pos[[sii[i]]]=c(insertion.cnt$pos[[sii[i]]], c(sis[i]:eis[i]))
        }
    }
    indel.cnt$pos=sapply(indel.cnt$pos, function(x) unique(x))
    indel.cnt$counts=as.vector(sapply(indel.cnt$pos, length))

    insertion.cnt$pos=sapply(insertion.cnt$pos, function(x) unique(x))
    insertion.cnt$counts=as.vector(sapply(insertion.cnt$pos, length))

    # this bit is modular and optional --------------------------------------------
    pam.indel.cnt = sapply(indel.cnt$pos, function(x) sum(x%in%pam.window.ind))
    pam.mm.cnt    = sapply(mm.cnt$pos,    function(x) sum(x%in%pam.window.ind))
   
    u1.indel.cnt = sapply(indel.cnt$pos, function(x) sum(x%in%u1.window.ind))
    u1.mm.cnt = sapply(mm.cnt$pos,    function(x) sum(x%in%u1.window.ind))

    u.indel.cnt = sapply(indel.cnt$pos, function(x) sum(x%in%u.window.ind))
    u.mm.cnt = sapply(mm.cnt$pos,    function(x) sum(x%in%u.window.ind))

    d1.indel.cnt = sapply(indel.cnt$pos, function(x) sum(x%in%d1.window.ind))
    d1.mm.cnt = sapply(mm.cnt$pos,    function(x) sum(x%in%d1.window.ind))

    d.indel.cnt = sapply(indel.cnt$pos, function(x) sum(x%in%d.window.ind))
    d.mm.cnt = sapply(mm.cnt$pos,    function(x) sum(x%in%d.window.ind))
    # ------------------------------------------------------------------------------

    mm.pos =sapply(mm.cnt$pos, paste, collapse=':')
    indel.pos=sapply(indel.cnt$pos, paste, collapse=':') 
    insertion.pos=sapply(insertion.cnt$pos, paste, collapse=':') 


        dftemp=data.frame(
        start.align=start(subject(pas)),
        end.align=end(subject(pas)),
        n.indel=indel.cnt$counts,
        n.insertion=insertion.cnt$counts,
        n.mismatch=nmismatch(pas),
        indel.pos=indel.pos,
        insertion.pos=insertion.pos,
        mm.pos=mm.pos, 
        pid=pid(pas),
        pam.indel.cnt=pam.indel.cnt,
        pam.mm.cnt= pam.mm.cnt, 
        u1.indel.cnt=u1.indel.cnt, 
        u1.mm.cnt=u1.mm.cnt, 
        u.indel.cnt=u.indel.cnt, 
        u.mm.cnt=u.mm.cnt, 
        d.indel.cnt=d.indel.cnt, 
        d.mm.cnt=d.mm.cnt, 
        d1.indel.cnt=d1.indel.cnt, 
        d1.mm.cnt=d1.mm.cnt,
        alignment=as.character(apas),
        stringsAsFactors=F)

    output= dftemp[match(observed.seqs, to.align.seqs),]
    oset[,grep(colnames(dftemp)[1], colnames(oset)):ncol(oset)]=output
    return(oset)
}
#giant.table=do.call('rbind', realigned.reads)
giant.table=data.frame(rbindlist(realigned.reads))
saveRDS(giant.table,  file = paste0(out.base.dir, 'processed/RData/giant.table.aligned2.RData'))

