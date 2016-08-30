library(rbamtools)
library(parallel)
library(lme4)

out.base.dir='/media/jbloom/d1/coupled_CRISPR/Experiments/082316/'

readID=list.files('/media/jbloom/d1/coupled_CRISPR/082316/fastq/')
# subsetting specific to this experiment
readID=unique(gsub('_R.*.fq.gz', '', readID))[1:12]

dir.create(paste(out.base.dir, 'processed/RData/', sep=''))



########################################################################################
load('/media/jbloom/d1/coupled_CRISPR/042916_Estops_Hiseq/ref/eStops_oligos.RData')
# 'oligos' contains information about each oligo

# how oligo names were constructed for reference contigs
oligo.name=apply(oligos,1,function(x) paste(x[c(1:3, 5:23, 25,26)], collapse=':'))

# flag oligos as coming from dubious ORFs-----------------------------------------------------------
load('/media/jbloom/d1/coupled_CRISPR/042916_Estops_Hiseq/ref/dubious_essential_genes.Rdata')
dubious.e.oligos=which(oligos$GENEID %in% dubious.essential)

oligos$dubious=FALSE
oligos$dubious[dubious.e.oligos]=TRUE
#---------------------------------------------------------------------------------------------------

# flag oligos that are actually causing synonymous changes -----------------------------------------
synolis=read.delim(file='/media/jbloom/d1/coupled_CRISPR/042916_Estops_Hiseq/ref/ESS_synonymous_controls.csv', header=T, sep='\t')
flagsyn=sort(as.vector(unlist(sapply(as.character(synolis$guide), function(x)     grep(x, oligos$guide)))))
oligos$flagsyn=FALSE
oligos$flagsyn[flagsyn]=TRUE
#---------------------------------------------------------------------------------------------------


# an experiment to hand check crispr induction for 10 transformed colonies 
hand.verified=read.delim('/media/jbloom/d1/coupled_CRISPR/042916_Estops_Hiseq/ref/hand_verified.txt', header=F, sep=':', stringsAsFactors=F)
checkme=lapply( hand.verified[,1], function(x) grep(x, oligos$stopPAM_oligos))
checkme=do.call('c', checkme)


# is there a memory leak somewhere ??? 2.6 GB of bams are blowing up memory usage
all.counts=list()
for( rr in readID) {
    print(rr)
    bam.file=paste(out.base.dir, 'processed/bam/', rr, '.bam', sep='')

    # use rbamtools to extract aligned results 
    reader<-bamReader(bam.file,idx=TRUE)
    contigs<-(refSeqDict(getHeaderText(reader)))@SN
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


save(all.counts, file=paste0(out.base.dir, '/processed/RData/all.counts.RData') )
load( paste0(out.base.dir, '/processed/RData/all.counts.RData') )
#load('/media/jbloom/d1/coupled_CRISPR/062916_Reworked/processed/RData/all.counts.RData')

#report non-zero counts
sapply(all.counts, function(x) sum(x$nAligns>0))

#total counts
sapply(all.counts, function(x) sum(x$nAligns))

# load all the raw data, intermediate structue, don't need to keep this
expt=list()
for(rr in readID) {  print(rr)
                     expt[[rr]]=readRDS(paste0(out.base.dir, 'processed/RData/', rr, '.RDS'))    }

names(expt)=readID
names(expt)=gsub('72-plate', '72plate', names(expt))
saveRDS(expt, file = paste0(out.base.dir, 'processed/RData/expt.RDS'))

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

swws.key=apply(expt.key, 1, function(x) paste( x[1],'_', x[6], sep=''))
wssw.key=apply(expt.key, 1, function(x) paste( x[3],'_', x[6], sep=''))

names(swws)=swws.key
names(wssw)=wssw.key

swws=swws[sort(names(swws))]
wssw=wssw[sort(names(wssw))]

c.expt=list()
for(i in 1:length(swws)) {
    c.expt[[names(swws)[i]]]=mapply('rbind', swws[[i]], wssw[[i]], SIMPLIFY=FALSE)
}
saveRDS(c.expt, file = paste0(out.base.dir, 'processed/RData/combined_expt.RDS'))


all.counts2 = sapply(c.expt, function(x) sapply(x, function(y) nrow(y) )
all.counts.perfect = sapply(c.expt, function(x) sapply(x, function(y) {
                                                           if(nrow(y)==0) {return(0) } 
                                                           else {
                                                              sum(y$cigar=='101M') }
                                                            } ))
all.barcode.diversity= sapply(c.expt, function(x) sapply(x, function(y) { 
                                                             if(nrow(y)==0) {return(0) } else {
                                                             length(unique(sort(y$barcode))) 
                                                            }} ))

t.test(log2(all.counts2[,'NMDm_96']+.5)~oligos$flagsyn)
t.test(log2(all.counts2[,'NMDm_96']+.5)~oligos$dubious)

t.test(log2(all.counts2[,'NMDm_72']+.5)~oligos$flagsyn)
t.test(log2(all.counts2[,'NMDm_72']+.5)~oligos$dubious)

t.test(log2(all.counts2[,'NMDm_72plate']+.5)~oligos$flagsyn)
t.test(log2(all.counts2[,'NMDm_72plate']+.5)~oligos$dubious)



t.test(log2(all.counts2[,'NMDm_0']+.5)~oligos$flagsyn)
t.test(log2(all.counts2[,'NMDm_0']+.5)~oligos$dubious)


t.test(log2(all.counts2[,'WT_96']+.5)~oligos$flagsyn)
t.test(log2(all.counts2[,'WT_96']+.5)~oligos$dubious)

t.test(log2(all.counts2[,'WT_72']+.5)~oligos$flagsyn)
t.test(log2(all.counts2[,'WT_72']+.5)~oligos$dubious)

t.test(log2(all.counts2[,'WT_72plate']+.5)~oligos$flagsyn)
t.test(log2(all.counts2[,'WT_72plate']+.5)~oligos$dubious)

t.test(log2(all.counts2[,'WT_0']+.5)~oligos$flagsyn)
t.test(log2(all.counts2[,'WT_0']+.5)~oligos$dubious)


# reads perfect
t.test(log2(all.counts.perfect[,'NMDm_96']+.5)~oligos$flagsyn)
t.test(log2(all.counts.perfect[,'NMDm_96']+.5)~oligos$dubious)

t.test(log2(all.counts.perfect[,'NMDm_72']+.5)~oligos$flagsyn)
t.test(log2(all.counts.perfect[,'NMDm_72']+.5)~oligos$dubious)

t.test(log2(all.counts.perfect[,'NMDm_72plate']+.5)~oligos$flagsyn)
t.test(log2(all.counts.perfect[,'NMDm_72plate']+.5)~oligos$dubious)


#--------------------------------------------------------------------------------------------------
t.test(log2(all.counts.perfect[,'NMDm_0']+.5)~oligos$flagsyn)
t.test(log2(all.counts.perfect[,'NMDm_0']+.5)~oligos$dubious)


t.test(log2(all.counts.perfect[,'WT_96']+.5)~oligos$flagsyn)
t.test(log2(all.counts.perfect[,'WT_96']+.5)~oligos$dubious)

t.test(log2(all.counts.perfect[,'WT_72']+.5)~oligos$flagsyn)
t.test(log2(all.counts.perfect[,'WT_72']+.5)~oligos$dubious)

t.test(log2(all.counts.perfect[,'WT_72plate']+.5)~oligos$flagsyn)
t.test(log2(all.counts.perfect[,'WT_72plate']+.5)~oligos$dubious)

t.test(log2(all.counts.perfect[,'WT_0']+.5)~oligos$flagsyn)
t.test(log2(all.counts.perfect[,'WT_0']+.5)~oligos$dubious)
#---------------------------------------------------------------------------------------------

# to do 
# change in oligos over time as a function of distance from end
wt.counts.perfect=all.counts.perfect[,c(7:10,12)]

wt.norm=t(t(wt.counts.perfect)/colSums(wt.counts.perfect))

par(mfrow=c(3,1))
plot(oligos$dist_from_CDS_end, wt.norm[,'WT_0'], ylim=c(0, 0.006))
plot(oligos$dist_from_CDS_end, wt.norm[,'WT_72'], ylim=c(0, 0.006))
plot(oligos$dist_from_CDS_end, wt.norm[,'WT_96'], ylim=c(0, 0.006))

identify(oligos$dist_from_CDS_end, wt.norm[,'WT_96'], oligos$G

wgs=split(data.frame(wt.norm), oligos$GENEID)
w96=sapply(wgs, function(x) x[,5])
w0=sapply(wgs, function(x) x[,1])

dist.from.end.list=split(oligos$dist_from_CDS_end, oligos$GENEID)


wfit=coefficients(lm(t(wt.norm)~t.var))[2,]
wfit.by.gene=split(wfit, oligos$GENEID)

w.b.e=mapply(cor, wfit.by.gene, dist.from.end.list)

dgenes=split(oligos$dubious, oligos$GENEID)
dgenes=sapply(dgenes, sum)>0



t.test(w.b.e~dgenes)


c1=mapply(cor, w96,dist.from.end.list, method='spearman')
c0=mapply(cor, w0,dist.from.end.list)

t.test(log2(all.counts.perfect[,'NMDm_96']+.5)~oligos$flagsyn)


t.test(log2(all.counts2[['106-295:131-296_hr_96']]$nAligns+.5)~oligos$flagsyn)
t.test(log2(all.counts2[['131-295:106-296_hr_72-plate']]$nAligns+.5)~oligos$dubious)
t.test(log2(all.counts2[['106-295:131-296_hr_72-plate']]$nAligns+.5)~oligos$dubious)



plot(oligos$dist_from_CDS_end/oligos$CDS_length , wt.norm[,5]/wt.norm[,1], ylim=c(0,2),col=ifelse(oligos$dubious,'red', 'black')    )

x2=oligos$dist_from_CDS_end/oligos$CDS_length 
y2=wt.norm[,5]/(wt.norm[,1])

finite.ind=is.finite(wt.norm[,5]/wt.norm[,1])
cor.test(c(oligos$dist_from_CDS_end/oligos$CDS_length)[finite.ind] , c(wt.norm[,5]/wt.norm[,1])[finite.ind])
cor.test(x2,y2)
plot(x2[finite.ind ],y2[finite.ind],col=ifelse(oligos$dubious[finite.ind],'red', 'black')    )
abline(lm( y2[finite.ind]~x2[finite.ind] ), col='blue')
       
par(mfrow=c(2,1))
plot(x2[finite.ind & !oligos$dubious & wt.counts.perfect[,1]>30  ],y2[finite.ind & !oligos$dubious & wt.counts.perfect[,1]>30], ylim=c(0,10) )#,col=ifelse(oligos$dubious[finite.ind],'red', 'black')    )
plot(x2[finite.ind & oligos$dubious & wt.counts.perfect[,1]>30 ],y2[finite.ind & oligos$dubious & wt.counts.perfect[,1]>30], ylim=c(0,10) )#,col=ifelse(oligos$dubious[finite.ind],'red', 'black')    )


cor.test(x2[finite.ind & !oligos$dubious & wt.counts.perfect[,1]>30  ],y2[finite.ind & !oligos$dubious & wt.counts.perfect[,1]>30])
       
       , wt.norm[,5]/wt.norm[,1], ylim=c(0,2),col=ifelse(oligos$dubious,'red', 'black')    )




t.test(x2~(y2>0.1 & wt.counts.perfect[,1]>30 & !oligos$dubious))


y3 =y2[finite.ind & !oligos$dubious & wt.counts.perfect[,1]>30 ]
gindex= oligos$GENEID[finite.ind & !oligos$dubious & wt.counts.perfect[,1]>30 ]

ygs=split(y3,gindex)
ygs.length=sapply(ygs, length)

plot(sapply(ygs, mean)[ygs.length>3], sapply(ygs,var)[ygs.length>3] , xlim=c(0,2), ylim=c(0,10))
identify(sapply(ygs, mean)[ygs.length>3], sapply(ygs,var)[ygs.length>3],  names(ygs)[ygs.length>3])

c.expt[[12]][oligos$GENEID=='YMR028W']

c.expt[[12]][oligos$GENEID=='YMR028W']

identify((sapply(ygs, mean)), (sapply(ygs,var)), 

         oligos[oligos$GENEID=='YOR151C',]
         sapply(c.expt[[7]][oligos$GENEID=='YOR151C'], nrow)
         sapply(c.expt[[12]][oligos$GENEID=='YOR151C'], nrow)


#by gene 



# 
# 
# mapply(c.expt, function(x) {lapply(x, function(y) {
#   if(length(y$barcode)>0 ) {
#     z=rle(sort(y$barcode))
#     z=data.frame(z$values, z$lengths)
#     z=z[order(z[,2], decreasing=T),]
#     return(z)
#     }
#     else {
#         return(NULL) # rbind(c('NNNNN', 0)))
#     }
# }) } )
# 
# 




cb.df=data.frame(t0=all.counts[[1]]$nAligns+all.counts[[7]]$nAligns,
                 t24=all.counts[[2]]$nAligns+all.counts[[8]]$nAligns,
                 t48=all.counts[[3]]$nAligns+all.counts[[9]]$nAligns,
                 t72=all.counts[[6]]$nAligns+all.counts[[11]]$nAligns,
                 t96=all.counts[[6]]$nAligns+all.counts[[12]]$nAligns)
cb.df=data.matrix(cb.df)
t.var=c(0,24,48,72,96)

cb.df.norm=t(t(cb.df)/colSums(cb.df))



tfit=lm(t(cb.df.norm)~t.var)
tnull=lm(t(cb.df.norm)~1)


 null= lm(t.tpm.matrix~covariates.OD)
   full= lm(t.tpm.matrix~covariates.OD+gdata[,hp])

rn=residuals(tnull)
fn=residuals(tfit) 
rss1=colSums(rn^2)
rss2=colSums(fn^2)
Fstat=(rss1-rss2)/(rss2/(5-1))
 pFstat=pf(Fstat, 1,5-1, lower.tail=F)
 cor.test(oligos$dist_from_CDS_end, -log10(pFstat))

plot(log10(oligos$dist_from_CDS_end), tfit$effects[2,])
identify(oligos$dist_from_CDS_end, tfit$effects[2,], oligos$GENEID)





#13 and 14 switch places with 15 and 16

plot(oligos$dist_from_CDS_end, all.counts[[12]]$nAligns)
plot(oligos$dist_from_CDS_end, all.counts[[1]]$nAligns)

plot(log2(all.counts[['plate_glu']]$nAligns+1), log2(all.counts[['plate_gal']]$nAligns+1), col =oligos$flagsyn+1 )
plot(log10(oligos$dist_from_CDS_end), log10(all.counts[[12]]$nAligns/(all.counts[[1]]$nAligns+.01)  ) )

cor.test(log10(oligos$dist_from_CDS_end), log10(all.counts[[12]]$nAligns/(all.counts[[1]]$nAligns+.01) +.5 ) )


t.test(log2(all.counts[['131-295:106-296_hr_96']]$nAligns+.5)~oligos$flagsyn)
t.test(log2(all.counts[['106-295:131-296_hr_96']]$nAligns+.5)~oligos$flagsyn)
t.test(log2(all.counts[['131-295:106-296_hr_72-plate']]$nAligns+.5)~oligos$dubious)
t.test(log2(all.counts[['106-295:131-296_hr_72-plate']]$nAligns+.5)~oligos$dubious)










for(i in 1:12) {
    print(readID[i])
    print(t.test(all.counts[[i]]$nAligns~oligos$flagsyn))
    print(t.test(all.counts[[i]]$nAligns~oligos$dubious))    
    #print(wilcox.test(all.counts[[i]]$nAligns~oligos$flagsyn))
    #print(wilcox.test(all.counts[[i]]$nAligns~oligos$dubious))
    #boxplot(log2(all.counts[[i]]$nAligns+1)~oligos$flagsyn)
    #boxplot(log2(all.counts[[i]]$nAligns+1)~oligos$dubious)
    #plot(jitter(as.numeric(oligos$dubious)), log2(all.counts[[i]]$nAligns+1))

}

(oligos[which(all.counts[[16]]$nAligns>100),c(-24,-27)])


plot(all.counts[[15]]$nAligns,all.counts[[16]]$nAligns)
plot(all.counts[[13]]$nAligns,all.counts[[14]]$nAligns)
plot(all.counts[[13]]$nAligns,all.counts[[5]]$nAligns)
plot(all.counts[[14]]$nAligns,all.counts[[6]]$nAligns)
plot(all.counts[[5]]$nAligns, all.counts[[6]]$nAligns)



plot(oligos$dist_from_CDS_end, all.counts[[13]]$nAligns-all.counts[[5]]$nAligns )
plot(oligos$dist_from_CDS_end, all.counts[['14']]$nAligns)
plot(oligos$dist_from_CDS_end, all.counts[['5']]$nAligns+all.counts[['6']]$nAligns)
plot(oligos$dist_from_CDS_end, all.counts[['6']]$nAligns)
plot(oligos$dist_from_CDS_end, all.counts[['13']]$nAligns/sum(all.counts[['13']]$nAligns)-all.counts[['5']]$nAligns/sum(all.counts[['5']]$nAligns ))
abline(h=0, col='red')

x=readRDS('/media/jbloom/d1/coupled_CRISPR/062916_Reworked/processed/RData/5.RDS')

z=sapply(x, function(y) y$barcode.type)
zz=unlist(z)
rle(sort(zz))

#load all.counts
#load('/media/jbloom/d1/coupled_CRISPR/042916_Estops_Hiseq/processed/RData/all.counts.RData')

gal.nmdpos=which(grepl('gal', names(all.counts)) & grepl('pos', names(all.counts)))
gal.nmdminus=which(grepl('gal', names(all.counts)) & grepl('min', names(all.counts)))

dex.nmdpos=which(grepl('dex', names(all.counts)) & grepl('pos', names(all.counts)))
dex.nmdminus=which(grepl('dex', names(all.counts)) & grepl('min', names(all.counts)))

gal.nmdpos.list=all.counts[c(dex.nmdpos[1], gal.nmdpos)]
gal.nmdminus.list=all.counts[c(dex.nmdminus[1], gal.nmdminus)]
dex.nmdpos.list=all.counts[dex.nmdpos]
dex.nmdminus.list=all.counts[dex.nmdminus]
plot(oligos$dist_from_CDS_end, all.counts[['plate_gal']]$nAligns)


t.var=c(0,4,8,24,72,96)

#double check 
all.combos=expand.grid(c(1:22), c(1:22))
for(n in 1:nrow(all.combos)) {
 print(   identical( rownames(all.counts[all.combos[n,1]]), rownames(all.counts[all.combos[n,2]])))
}

dm=sapply(dex.nmdminus.list, function(x) x$nAligns)
dp=sapply(dex.nmdpos.list, function(x) x$nAligns)
gm=sapply(gal.nmdminus.list, function(x) x$nAligns)
gp=sapply(gal.nmdpos.list, function(x) x$nAligns)





#t(t(gp[order(gp[,6], decreasing=T)[1:20],])/colSums(gp))
t(t(gm[order(gm[,6], decreasing=T)[1:20],])/colSums(gm))

t(t(dp[order(dp[,6], decreasing=T)[1:20],])/colSums(dp))
t(t(dm[order(dm[,6], decreasing=T)[1:20],])/colSums(dm))

#GP
gp[order(gp[,6], decreasing=T)[1:20],6]
"YDR160W" "YKR004C" "E  YAL025C (near middle)"  "E YBR256C (near end)" 'E YAL043C (near beginning)'   "E (near end) YLR259C"


x=readRDS(file=paste0(out.base.dir, 'processed/RData/', readID[1], '.RDS'))

1757





#build RDS structures by experiment-----------------------------------------------------------------
buildRDS.expt=function(rin,fout) {
  expt=list()
    for ( rr in rin) {
        print(rr)
        expt[[rr]]=readRDS(paste0(out.base.dir, 'processed/RData/', rr, '.RDS'))
    }
    saveRDS(expt, file=fout)
}
rin=readID[c(dex.nmdpos[1], gal.nmdpos)]
fout=paste0(out.base.dir, 'processed/RData/', 'galpos.RDS')
buildRDS.expt(rin, fout)

rin=readID[c(dex.nmdminus[1], gal.nmdminus)]
fout=paste0(out.base.dir, 'processed/RData/', 'galminus.RDS')
buildRDS.expt(rin, fout)

rin=readID[dex.nmdpos]
fout=paste0(out.base.dir, 'processed/RData/', 'dexpos.RDS')
buildRDS.expt(rin, fout)

rin=readID[dex.nmdminus]
fout=paste0(out.base.dir, 'processed/RData/', 'dexminus.RDS')
buildRDS.expt(rin, fout)
#----------------------------------------------------------------------------------------------------

galpos.RDS=paste0(out.base.dir, 'processed/RData/', 'galpos.RDS')
galminus.RDS=paste0(out.base.dir, 'processed/RData/', 'galminus.RDS')
dexpos.RDS=paste0(out.base.dir, 'processed/RData/', 'dexpos.RDS')
dexminus.RDS=paste0(out.base.dir, 'processed/RData/', 'dexminus.RDS')



expt=readRDS(galpos.RDS)

c.b.bc=lapply(expt, function(x){ lapply(x, function(y) {
    if(length(y$barcode)>0 ) {
    z=rle(sort(y$barcode))
    z=data.frame(z$values, z$lengths)
    z=z[order(z[,2], decreasing=T),]
    return(z)
    }
    else {
        return(NULL) # rbind(c('NNNNN', 0)))
    }
}) } )

cb=lapply(c.b.bc, function(x) sapply(x, function(y) { if(is.null(y)) {return(0)} else {nrow(y)} }))
cb.df=do.call('cbind', cb)                   
                   
                                 rle


extracted.reads=lapply(expt, function(x) x[checkme])
outtable=list()
counter=1
for(n in names(extracted.reads[[1]]) ) {
    print(n)
    outtable[[n]]= ldply(extracted.reads, function(x) x[[n]])
    write.table(outtable[[n]], file=paste0(out.base.dir, 'processed/RData/galpos_', counter, '.txt'), sep='\t', quote=F)
    counter=counter+1
    if(counter==3) {counter=5}
}

expt=readRDS(galminus.RDS)
extracted.reads=lapply(expt, function(x) x[checkme])
outtable=list()
counter=1
for(n in names(extracted.reads[[1]]) ) {
    print(n)
    outtable[[n]]= ldply(extracted.reads, function(x) x[[n]])
    write.table(outtable[[n]], file=paste0(out.base.dir, 'processed/RData/galminus_', counter, '.txt'), sep='\t', quote=F)
    counter=counter+1
    if(counter==3) {counter=5}
}

expt=readRDS(dexpos.RDS)
extracted.reads=lapply(expt, function(x) x[checkme])
outtable=list()
counter=1
for(n in names(extracted.reads[[1]]) ) {
    print(n)
    outtable[[n]]= ldply(extracted.reads, function(x) x[[n]])
    write.table(outtable[[n]], file=paste0(out.base.dir, 'processed/RData/dexpos_', counter, '.txt'), sep='\t', quote=F)
    counter=counter+1
    if(counter==3) {counter=5}
}

expt=readRDS(dexminus.RDS)
extracted.reads=lapply(expt, function(x) x[checkme])
outtable=list()
counter=1
for(n in names(extracted.reads[[1]]) ) {
    print(n)
    outtable[[n]]= ldply(extracted.reads, function(x) x[[n]])
    write.table(outtable[[n]], file=paste0(out.base.dir, 'processed/RData/dexminus_', counter, '.txt'), sep='\t', quote=F)
    counter=counter+1
    if(counter==3) {counter=5}
}






# extract only the perfect matching reads for each experiment --------------------------------------
gp.good=gp
rin=readID[c(dex.nmdpos[1], gal.nmdpos)]
for ( rr in rin) {
    print(rr)
    br=readRDS(paste0(out.base.dir, 'processed/RData/', rr, '.RDS'))
    gp.good[,rr]=sapply(br, function(x) sum(x$cigar=='101M'))
}
readID[c(dex.nmdminus[1], gal.nmdminus)]
readID[dex.nmdpos]
readID[dex.nmdminus]
#---------------------------------------------------------------------------------------------------


#correlation between all counts and counts coming from perfect matches increases at later time points

#calc poisson regression
calc_pois_regression=function(counts, t.var) {
    tot.counts=colSums(counts)
    #bad.exp=c(5:6) #which(tot.counts<2e6)
    #counts[,bad.exp]=NA
    nm=mcMap(function(i){
                    glm(counts[i,]~offset(log(tot.counts)), family='poisson')
                },   1:nrow(counts),mc.cores=30)
    fm=mcMap(function(i){
                    glm(counts[i,]~t.var+offset(log(tot.counts)), family='poisson')
                },   1:nrow(counts),mc.cores=30)
 
    am=simplify2array(mcMap(function(i){
               anova(fm[[i]],nm[[i]], test='Chisq')$Deviance[2] 
                }, 1:nrow(counts), mc.cores=30))
    return(list(nm=nm, fm=fm,am=am))
}


cb.df=data.frame(t0=all.counts[[1]]$nAligns+all.counts[[7]]$nAligns,
                 t24=all.counts[[2]]$nAligns+all.counts[[8]]$nAligns,
                 t48=all.counts[[3]]$nAligns+all.counts[[9]]$nAligns,
                 t72=all.counts[[6]]$nAligns+all.counts[[11]]$nAligns,
                 t96=all.counts[[6]]$nAligns+all.counts[[12]]$nAligns)
cb.df=data.matrix(cb.df)
t.var=c(0,24,48,72,96)




cb.fit=calc_pois_regression(cb.df, t.var)
cb.c=sapply(cb.fit$fm, function(x) coef(x))
efit=exp(cb.c[2,])
efit2=scale(efit)
#oligos$dubious
m0=lmer(cb.c[2,]~1+oligos$flagsyn+oligos$dist_from_CDS_end+(1|as.factor(oligos$GENEID)))

m1=lmer(cb.c[2,]~1+(1|as.factor(oligos$GENEID)))


dm.fit=calc_pois_regression(dm, t.var)
dp.fit=calc_pois_regression(dp, t.var)
gm.fit=calc_pois_regression(gm, t.var)

gp.fit=calc_pois_regression(gp, t.var)
gp.c=sapply(gp.fit$fm, function(x) coef(x))
efit=exp(gp.c[2,])
efit2=scale(efit)
m0=lmer(efit2~1+(1|as.factor(oligos$GENEID)))
m1=lmer(efit2~1+oligos$flagsyn+(1|as.factor(oligos$GENEID)))


gpg.fit=calc_pois_regression(gp.good, t.var)
gpg.c=sapply(gpg.fit$fm, function(x) coef(x))
efitg=exp(gpg.c[2,])
efit2g=scale(efitg)
m0g=lmer(efit2g~1+(1|as.factor(oligos$GENEID)))
m1g=lmer(efit2g~1+oligos$dubious+(1|as.factor(oligos$GENEID)))

gT=gp.good[,1]/sum(gp.good[,1])
gE=gp.good[,6]/sum(gp.good[,6])

plot(oligos$dist_from_CDS_end/oligos$CDS_length, gE/gT, ylim=c(0,.1))

wilcox.test(gE/gT~oligos$dubious)


gT=gp[,1]/sum(gp[,1])
gE=gp[,6]/sum(gp[,6])

dropouts=which(gE/gT==0)
dropped= 1:10971 %in% dropouts
len=oligos$dist_from_CDS_end/oligos$CDS_length

anova(m0,m1)
o2=oligos
o2$efit=efit
o3=split(o2, o2$GENEID)
plot(o3[[2]]$dist_from_CDS_end, o3[[2]]$efit)

m1=lm(efit~as.factor(oligos$GENEID)+oligos$dist_from_CDS_end+as.factor(oligos$guideStrand)+as.factor(oligos$case) + as.factor(oligos$codingStrand))
anova(m1)
m0=lmer(efit~1 +(1|1))

efit2=scale(efit)
m0=lmer(efit2~1+(1|as.factor(rep(1,10971))) #(1|as.factor(oligos$GENEID)))

rstat=VarCorr(m1)[[1]][1,1]
kstat=rep(0,1000)
for(k in 376:1000){
    print(k)
    m1=lmer(sample(efit2)~1+(1|as.factor(oligos$GENEID)))
    kstat[k]=VarCorr(m1)[[1]][1,1]
}

anova(m0,m1)
sp=stack(p)  
   reffMod = lmer(z ~ 1 + (1|as.factor(id)))
    VarCorr(reffMod)
     Vcomp = as.numeric((attributes(summary(reffMod)))$REmat[,'Variance'])
     VarG = Vcomp[1]
     VarE = Vcomp[2]
     H2 = VarG/(VarE+VarG)
regress(z~1,~as.factor(id))

#m1=lm(efit~as.factor(oligos$guideStrand)+as.factor(oligos$case) + as.factor(oligos$codingStrand))
#anova(m1)

m1.p=lm(sample(efit)~as.factor(oligos$GENEID)+oligos$dist_from_CDS_end)
anova(m1.p)

aDev=mcmapply(function(fm1, nm1) {
        anova(fm1,nm1, test='Chisq')$Deviance[2]
}, nm, fm, mc.cores=11)

oligos$GENEID[which(efit2>1)]

ii

    nm[[i]]=glm(counts[i,]~offset(log(tot.counts)), family='poisson')
    fm[[i]]=glm(counts[i,]~t.var+offset(log(tot.counts)), family='poisson' )





source("https://bioconductor.org/biocLite.R")
library("AnnotationDbi")
library("org.Sc.sgd.db")





dm=sapply(dex.nmdminus.list, function(x) x$nAligns)
dm.norm=dm/colSums(dm)
dp=sapply(dex.nmdpos.list, function(x) x$nAligns)
dp.norm=dp/colSums(dp)

m0=lm(t(dp.norm)~1)
m1=lm(t(dp.norm)~t.var)
hist(coefficients(m1)[2,], breaks=100000, xlim=c(-.00002, .00002), main='dextrose nmd +')
rn=residuals(m0)
fn=residuals(m1) 
rss1=colSums(rn^2)
rss2=colSums(fn^2)
Fstat=((rss1-rss2)/(2-1))/(rss2/(6-2))
p.fit=pf(Fstat, 1, 6-1, lower.tail=F)
x11()
hist(p.fit, breaks=1000, main='dextrose nmd +')

gp=sapply(gal.nmdpos.list, function(x) x$nAligns)
rownames(gp)=contigs
library(NBPSeq)


mm=model.matrix(lm(gp[1,]~t.var))
beta0=c(NA,0)
model.test=nb.glm.test(gp, mm, beta0, normalization.method=NULL, subset=1:1000)
norm.factors=estimate.norm.factors(gp)
norm.factors=estimate.norm.factors(gp, method=NULL)

gp.nb=prepare.nb.data(gp)i
dest=estimate.disp(obj=gp.nb)

test.coefficient
library(MASS)
plot(gp[1,])
plot(gp[1,])

tot.counts=colSums(gp)
nm=list()
fm=list()
mcomp=list()
#why is 5012 failing
#7362 (all 0)
#10344 (all 0)
for(i in 1:10971) {
    print(i)
    #nm[[i]]=glm.nb(gp[i,]~offset(log(tot.counts)))
    #fm[[i]]=glm.nb(gp[i,]~t.var+offset(log(tot.counts)))
    nm[[i]]=glm(gp[i,]~offset(log(tot.counts)), family='poisson' )
    fm[[i]]=glm(gp[i,]~t.var+offset(log(tot.counts)), family='poisson' )
}

sapply(mcomp, function(x) x$'Pr(Chi)'[2])
hist(do.call('c', sapply(fm, function(x) coef(x)[2])), breaks=1000, xlim=c(-.3,.1))
nbfit=sapply(fm, function(x) {
             cout=coef(x)[2] 
             if(is.null(cout))
             cout=0
             return(cout)}
             )
plot(log10(y+1), nbfit, ylim=c(-.2,.1))

exp(predict(glm.nb(gp[202,]~t.var+offset(log(tot.counts)))))
exp(predict(glm(gp[202,]~t.var+offset(log(tot.counts)), family='poisson' )))



anova(nm[[5]], fm[[5]])

#glmNullWithReplicates <- glmer.nb(thisRNA ~ log(thisDNA) + (1|thisBarcode) + (1|thisReplicate), data=glmMatrixWithReplicates)
#glmer.nb(gp[101,]~t.var+log(tot.counts) + (1|gp[101,]))

tfit2=lm(gp[100,]~t.var+offset(log(tot.counts)))
plot(gp[1,], predict(tfit2))
plot(gp[1,], exp(predict(tfit)))


gp.norm=gp/colSums(gp)
apply(gp, 2, function(x) sum(x==0))
gm=sapply(gal.nmdminus.list, function(x) x$nAligns)
gm.norm=gm/colSums(gm)
apply(gm, 2, function(x) sum(x==0))

m0=lm(log(t(gp.norm)+1)~1)
m1=lm(log(t(gp.norm)+1)~t.var)
gpc=coefficients(m1)[2,]

m0=lm(t(gm.norm)~1)
m1=lm(t(gm.norm)~t.var)
gmc=coefficients(m1)[2,]


x11()
hist(coefficients(m1)[2,], breaks=100000, xlim=c(-.00002, .00002),main='gal nmd +')
rn=residuals(m0)
fn=residuals(m1) 
rss1=colSums(rn^2)
rss2=colSums(fn^2)
Fstat=((rss1-rss2)/(2-1))/(rss2/(6-2))
p.fit=pf(Fstat, 1, 6-1, lower.tail=F)
x11()
hist(p.fit, breaks=1000, main='gal nmd +')  




m0=lm(t(gm.norm)~1)
m1=lm(t(gm.norm)~t.var)
x11()
hist(coefficients(m1)[2,], breaks=100000, xlim=c(-.00002, .00002),main='gal nmd -')
rn=residuals(m0)
fn=residuals(m1) 
rss1=colSums(rn^2)
rss2=colSums(fn^2)
Fstat=((rss1-rss2)/(2-1))/(rss2/(6-2))
p.fit=pf(Fstat, 1, 6-1, lower.tail=F)
x11()
hist(p.fit, breaks=1000, main='gal nmd -')  


m0=lm(t(dm.norm)~1)
m1=lm(t(dm.norm)~t.var)
x11()
hist(coefficients(m1)[2,], breaks=100000, xlim=c(-.00002, .00002),main='dex nmd -')
rn=residuals(m0)
fn=residuals(m1) 
rss1=colSums(rn^2)
rss2=colSums(fn^2)
Fstat=((rss1-rss2)/(2-1))/(rss2/(6-2))
p.fit=pf(Fstat, 1, 6-1, lower.tail=F)
x11()
hist(p.fit, breaks=1000, main='dex nmd -')  


x=coefficients(m1)[2,]

y=oligos$dist_from_CDS_end


plot(x,log10(y), xlim=c(-.00002, .00002), col='#00000011', pch=20)

syn_controls=read.delim(file='/media/jbloom/d1/coupled_CRISPR/042916_Estops_Hiseq/ref/ESS_synonymous_controls.csv',
                        sep='\t', header=T, stringsAsFactors=F)
syc.ind=match(syn_controls$Index, oligos$Index)
points(x[syc.ind],log10(y)[syc.ind], xlim=c(-.00002, .00002), col='#ff000066', pch=20)

which(x>0)

load('/media/jbloom/d1/coupled_CRISPR/042916_Estops_Hiseq/ref/dubious_essential_genes.Rdata')
dubious.e.oligos=which(oligos$GENEID %in% dubious.essential)



which(nbfit>0) %in% dubious.e.oligos


match(dubious.essential, oligos$GENEID) 
o

apply(dm,2, function(x) (x>0))

gal.pos.dropouts=which(gp[,6]==0)


t0=readRDS('/media/jbloom/d1/coupled_CRISPR/042916_Estops_Hiseq/processed/RData/1-dex_0_nmdpos.RDS')
ubc0=sapply(t0, function(y) 
           tryCatch( {(unique(sort(y$barcode))) }, error=function(e) return(0)                  
                    ) )

ic0=sapply(ubc0, function(x) sum(sapply(x, function(y) is.character(y))))
#12353

t1=readRDS('/media/jbloom/d1/coupled_CRISPR/042916_Estops_Hiseq/processed/RData/7-gal_8_nmdpos.RDS')
ubc1=sapply(t1, function(y) 
           tryCatch( {(unique(sort(y$barcode))) }, error=function(e) return(0)                  
                    ) )

ic1=sapply(ubc1, function(x) sum(sapply(x, function(y) is.character(y))))


t2=readRDS('/media/jbloom/d1/coupled_CRISPR/042916_Estops_Hiseq/processed/RData/19-gal_96_nmdpos.RDS')
ubc2=sapply(t2, function(y) 
           tryCatch( {(unique(sort(y$barcode))) }, error=function(e) return(0)                  
                    ) )

ic2=sapply(ubc2, function(x) sum(sapply(x, function(y) is.character(y))))
#12353






ubc=lapply(t2, function(y) 
           tryCatch( {rle((sort(y$barcode))) }, error=function(e) return(0)                  
                    ) )
glmNullWithReplicates <- glmer.nb(thisRNA ~ log(thisDNA) + (1|thisBarcode) + (1|thisReplicate), data=glmMatrixWithReplicates)

glmH1WithReplicates <- glmer.nb(thisRNA ~ log(thisDNA) + thisBYRM + (1|thisBarcode) + (1|thisReplicate), data=glmMatrixWithReplicates)

