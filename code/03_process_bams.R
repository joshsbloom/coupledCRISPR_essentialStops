library(rbamtools)
library(parallel)
library(lme4)

out.base.dir='/media/jbloom/d1/coupled_CRISPR/Experiments/082316/'

readID=list.files(paste0(out.base.dir, 'fastq/'))
# subsetting specific to this experiment
readID=unique(gsub('_R.*.fq.gz', '', readID))[1:12]

# output here 
dir.create(paste(out.base.dir, 'processed/RData/', sep=''))


# Load various annotation files ----------------------------------------------------------------------------------------
load('/media/jbloom/d1/coupled_CRISPR/Reference/eStops_oligos.RData')
# 'oligos' contains information about each oligo

# how oligo names were constructed for reference contigs
oligo.name=apply(oligos,1,function(x) paste(x[c(1:3, 5:23, 25,26)], collapse=':'))

# flag oligos as coming from dubious ORFs-----------------------------------------------------------
load('/media/jbloom/d1/coupled_CRISPR/Reference/dubious_essential_genes.Rdata')
dubious.e.oligos=which(oligos$GENEID %in% dubious.essential)

oligos$dubious=FALSE
oligos$dubious[dubious.e.oligos]=TRUE
#---------------------------------------------------------------------------------------------------

# flag oligos that are actually causing synonymous changes -----------------------------------------
synolis=read.delim(file='/media/jbloom/d1/coupled_CRISPR/Reference/ESS_synonymous_controls.csv', header=T, sep='\t')
flagsyn=sort(as.vector(unlist(sapply(as.character(synolis$guide), function(x)     grep(x, oligos$guide)))))
oligos$flagsyn=FALSE
oligos$flagsyn[flagsyn]=TRUE
#---------------------------------------------------------------------------------------------------


# an experiment to hand check crispr induction for 10 transformed colonies 
#hand.verified=read.delim('/media/jbloom/d1/coupled_CRISPR/Reference/hand_verified.txt', header=F, sep=':', stringsAsFactors=F)
#checkme=lapply( hand.verified[,1], function(x) grep(x, oligos$stopPAM_oligos))
#checkme=do.call('c', checkme)

#save(oligos, file=paste0(paste(out.base.dir, 'processed/RData/', sep=''), 'oligoAnnotations.RData'))
# ---------------------------------------------------------------------------------------------------------------------
# load oligo annotations 
load(paste0(paste(out.base.dir, 'processed/RData/', sep=''), 'oligoAnnotations.RData'))


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
    swws=mclapply(expt, function(x){ lapply(x, function(y) { return(y[y$barcode.type=='swws',])})}, mc.cores=66)
    wssw=mclapply(expt, function(x){ lapply(x, function(y) { return(y[y$barcode.type=='wssw',])})}, mc.cores=66)


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
    c.expt[[names(swws)[i]]]=mapply('rbind', swws[[i]], wssw[[i]], SIMPLIFY=FALSE)
}
#saveRDS(c.expt, file = paste0(out.base.dir, 'processed/RData/combined_expt.RDS'))

# read in list structure
c.expt=readRDS(file = paste0(out.base.dir, 'processed/RData/combined_expt.RDS'))

# all.counts2 (for combined experiment, total number of counts per expected oligo)
all.counts2 = sapply(c.expt, function(x) sapply(x, function(y) nrow(y) ))

# all.counts.perfect (for combined experiment, total number of counts per expected oligo where repair template perfectly matches expected sequence)
all.counts.perfect = sapply(c.expt, function(x) sapply(x, function(y) {
                                                           if(nrow(y)==0) {return(0) } 
                                                           else {
                                                              sum(y$cigar=='101M') }
                                                            } ))
all.barcode.diversity= sapply(c.expt, function(x) sapply(x, function(y) { 
                                                             if(nrow(y)==0) {return(0) } else {
                                                             length(unique(sort(y$barcode))) 
                                                            }} ))
count.list=list(all.counts2=all.counts2, all.counts.perfect=all.counts.perfect, all.barcode.diversity=all.barcode.diversity)
save(count.list, file = paste0(out.base.dir, 'processed/RData/count_tables.RData'))

#Exploratory analysis 

# t-tests (differences in counts per time point given likely synonymous or dubious status)
t.test.flagsyn.all=apply(log2(all.counts2+.5), 2, function(x) t.test(x~oligos$flagsyn))
t.test.flagsyn.perfect=apply(log2(all.counts.perfect+.5), 2, function(x) t.test(x~oligos$flagsyn))

t.test.dubious.all=apply(log2(all.counts2+.5), 2, function(x) t.test(x~oligos$dubious))
t.test.dubious.perfect=apply(log2(all.counts.perfect+.5), 2, function(x) t.test(x~oligos$dubious))

# extract all perfect counts from WT experiment (excluding 
wt.counts.perfect=all.counts.perfect[,c(7:10,12)]
wt.norm=t(t(wt.counts.perfect)/colSums(wt.counts.perfect))

par(mfrow=c(5,1))
plot(oligos$dist_from_CDS_end, wt.norm[,'WT_0'], ylim=c(0, 0.006))
plot(oligos$dist_from_CDS_end, wt.norm[,'WT_24'], ylim=c(0, 0.006))
plot(oligos$dist_from_CDS_end, wt.norm[,'WT_48'], ylim=c(0, 0.006))
plot(oligos$dist_from_CDS_end, wt.norm[,'WT_72'], ylim=c(0, 0.006))
plot(oligos$dist_from_CDS_end, wt.norm[,'WT_96'], ylim=c(0, 0.006))

#identify(oligos$dist_from_CDS_end, wt.norm[,'WT_96'], oligos$GENEID)

t.var=c(0,24,48,72,96)
plot(t.var, wt.norm[1,]*1e6, type='n', ylim=c(0,400),ylab='normalized counts * 1e6')
for(i in 1:10971) {  
    points(t.var, wt.norm[i,]*1e6, type='l', col='#00000011')  
}

# exploration (filter on total counts at t0)
par(mfrow=c(2,1))
plot(t.var, wt.norm[1,]*1e6, type='n', ylim=c(0,500),ylab='normalized counts * 1e6', xlab='time', main='not dubious')
t0count_gt30_notDubious = which(wt.counts.perfect[,1]>30 & !oligos$dubious)
for(i in t0count_gt30_notDubious) {  
    points(t.var, wt.norm[i,]*1e6, type='l', col='#00000011')  
}

plot(t.var, wt.norm[1,]*1e6, type='n', ylim=c(0,500),ylab='normalized counts * 1e6', xlab='time', main='dubious')
t0count_gt30_Dubious = which(wt.counts.perfect[,1]>30 & oligos$dubious)
for(i in t0count_gt30_Dubious) {  
    points(t.var, wt.norm[i,]*1e6, type='l', col='#00000011')  
}



# junk exploratory code (needs cleanup) --------------------------------------------------------------------------------------------------------
# plot change in normalized oligo count by length
plot(oligos$dist_from_CDS_end/oligos$CDS_length , wt.norm[,5]/wt.norm[,1], ylim=c(0,2),col=ifelse(oligos$dubious,'red', 'black')    )

# separate variables to keep track of division by zero
x2=oligos$dist_from_CDS_end/oligos$CDS_length 
y2=wt.norm[,5]/(wt.norm[,1])

finite.ind=is.finite(wt.norm[,5]/wt.norm[,1])
cor.test(c(oligos$dist_from_CDS_end/oligos$CDS_length)[finite.ind] , c(wt.norm[,5]/wt.norm[,1])[finite.ind])
cor.test(x2,y2)
plot(x2[finite.ind ],y2[finite.ind],col=ifelse(oligos$dubious[finite.ind],'red', 'black') , ylim=c(0,10)  )
abline(lm( y2[finite.ind]~x2[finite.ind] ), col='blue')
       
par(mfrow=c(2,1))
plot(x2[finite.ind & !oligos$dubious & wt.counts.perfect[,1]>30  ],y2[finite.ind & !oligos$dubious & wt.counts.perfect[,1]>30], ylim=c(0,10), 
     xlab='stop distance/gene length',
     ylab='normalized counts',
     main='not dubious' )#,col=ifelse(oligos$dubious[finite.ind],'red', 'black')    )
plot(x2[finite.ind & oligos$dubious & wt.counts.perfect[,1]>30 ],y2[finite.ind & oligos$dubious & wt.counts.perfect[,1]>30], ylim=c(0,10),
     xlab='stop distance/gene length',
     ylab='normalized counts',
     main='dubious' )#,col=ifelse(oligos$dubious[finite.ind],'red', 'black')    )

cor.test(x2[finite.ind & !oligos$dubious & wt.counts.perfect[,1]>30  ],y2[finite.ind & !oligos$dubious & wt.counts.perfect[,1]>30])
t.test(x2~(y2>0.1 & wt.counts.perfect[,1]>30 & !oligos$dubious))

# look at means and variances for oligo change ratio for each gene (excluding dubious orfs)
y3 =y2[finite.ind & !oligos$dubious & wt.counts.perfect[,1]>30 ]
gindex= oligos$GENEID[finite.ind & !oligos$dubious & wt.counts.perfect[,1]>30 ]

ygs=split(y3,gindex)
ygs.length=sapply(ygs, length)

plot(sapply(ygs, mean)[ygs.length>3], sapply(ygs,var)[ygs.length>3])
identify(sapply(ygs, mean)[ygs.length>3], sapply(ygs,var)[ygs.length>3],  names(ygs)[ygs.length>3])
c.expt[[12]][oligos$GENEID=='YMR028W']
c.expt[[12]][oligos$GENEID=='YMR028W']
#---------------------------------------------------------------------------------------------------------------------------------------------
