library(rbamtools)
library(parallel)
library(lme4)
library(stringdist)
library(mgcv)
library(bio3d)
library(msm)
library(foreach)
library(doParallel)
library(topGO)
library(org.Sc.sgd.db)
library(fdrtool)
library(sisus)
library(Biostrings)
library(stringr)
library(seqinr)
library(qvalue)
library(intervals)
library(data.table)
library(yeastExpData)
library(gdata)
library(lmerTest)
library(gamm4)
library(afex)
library(optimx)
library(car)

library("BSgenome.Scerevisiae.UCSC.sacCer3")
sacCer3=BSgenome.Scerevisiae.UCSC.sacCer3
library("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene")
txdb=TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
library(VariantAnnotation)
library(WriteXLS)
X11.options(type='cairo')

library(foreach)
library(doMC)
registerDoMC(cores=50)

x <- org.Sc.sgdALIAS
mapped_probes <- mappedkeys(x)
gene2alias <- sapply(as.list(x[mapped_probes]), paste, collapse=' ')


t.var=c(0,24,48,72,96)
out.base.dir='/media/jbloom/d1/coupled_CRISPR/Experiments/082316/'

# output here 
#dir.create(paste(out.base.dir, 'processed/RData/', sep=''))
# load oligos and annotation 
source('/media/jbloom/d1/coupled_CRISPR/code/accessory_functions.R')

# START HERE !!!!! 
seq.tables=readRDS(file = paste0(out.base.dir, 'processed/RData/seq.tables.RData'))
#giant.table=readRDS(file = paste0(out.base.dir, 'processed/RData/giant.table.RData'))
giant.table=readRDS(file = paste0(out.base.dir, 'processed/RData/giant.table.aligned2.RData'))

# minimum number of reads at t0
cutoff_t0=20
giant.table=addFilters(giant.table,cutoff_t0)

#Diagnostics Plots  --------------------------------------------------------------------------------------------------------------------
    par(mfrow=c(1,2))
    hist(giant.table$start.align, breaks=1000, xlim=c(0,20))
    hist(giant.table$end.align, breaks=1000, xlim=c(80,101))

    par(mfrow=c(2,1))
    hist(log10(giant.table$WT_0[giant.table$WT_0>0]), breaks=100, xlim=c(0,4),  main='WT t=0',
         xlab='log10 counts per barcode') 
    legend('topright', paste0(round(sum(giant.table$WT_0[giant.table$WT_0>0])/1e6,2),' x 1e6 total reads \n',
                       paste0(round(sum(giant.table$WT_0[giant.table$WT_0>cutoff_t0])/1e6,2),' x 1e6 total reads retained')))     
    abline(v=log10(cutoff_t0))
    hist(log10(giant.table$NMDm_0[giant.table$NMDm_0>0]), breaks=100, xlim=c(0,4), main='NMD t=0',
         xlab='log10 counts per barcode') 
    legend('topright', paste0(round(sum(giant.table$NMDm_0[giant.table$NMDm_0>0])/1e6,2),' x 1e6 total reads \n',
                       paste0(round(sum(giant.table$NMDm_0[giant.table$NMDm_0>cutoff_t0])/1e6,2),' x 1e6 total reads retained')))     
    abline(v=log10(cutoff_t0))

    #plot(ecdf(table(giant.table$WT_0[giant.table$WT_0>0])), xlim=c(0,100))

    expt.columns=grep('^WT|^NMD', names(giant.table))
    expt.counts=colSums(giant.table[,expt.columns])
    expt.counts.red=colSums(giant.table[giant.table$WT_0> cutoff_t0 | giant.table$NMDm_0> cutoff_t0,expt.columns])
    expt.counts.red2=colSums(giant.table[(giant.table$WT_0> cutoff_t0 | giant.table$NMDm_0> cutoff_t0) & giant.table$d1.filter  ,expt.columns])

    barplot(rbind(expt.counts.red2, expt.counts.red-expt.counts.red2, expt.counts-expt.counts.red),
            ylab='total reads', xlab='sample', main='reads per sample', col=c('black', 'grey', 'white'))
    legend('topright', c( paste0(round(sum(expt.counts)/1e6,1), ' x 1e6 total reads'),
                          paste0(round(sum(expt.counts.red)/1e6,1), ' x 1e6 total reads\n (reads per barcode filter)')  ,
                          paste0(round(sum(expt.counts.red2)/1e6,1), ' x 1e6 total reads\n (reads per barcode and perfect match repair)')),
                fill=c('black', 'grey', 'white'),
                   cex=1)  
#-------------------------------------------------------------------------------------------------------------------------------------

# Convert counts to RPKM (for plotting only)
giant.table.RPKM=10^9*t(t(giant.table[,expt.columns])/(101*expt.counts))
colnames(giant.table.RPKM)=paste0(colnames(giant.table.RPKM), '_rpkm')
giant.table=cbind(giant.table, giant.table.RPKM)
# identify barcodes that were identified in WT and NMD experiment 
giant.table$matched.barcode=giant.table$WT_0>cutoff_t0 & giant.table$NMDm_0>cutoff_t0

#215449 #R> sum(giant.table$cigar=='101M') #[1] 348288

# filter and rearrange table  
    giant.table.WT  = giant.table[giant.table$WT_0>cutoff_t0,]
    giant.table.WT  = giant.table.WT[,-which(grepl('^NMD', colnames(giant.table.WT)))]
    giant.table.NMD = giant.table[giant.table$NMDm_0>cutoff_t0,]
    giant.table.NMD  = giant.table.NMD[,-which(grepl('^WT', colnames(giant.table.NMD)))]

# apply filter for synthesis errors in repair templates 
    gt.WT.p   = giant.table.WT[giant.table.WT$u1_d_d1.filter,] #giant.table.WT$cigar=='101M',]
    gt.NMD.p  = giant.table.NMD[giant.table.NMD$u1_d_d1.filter,] #giant.table.NMD$cigar=='101M',]

# for example, find oligos with indel near cutting site 
    #gt.WT.del=(giant.table.WT[grep('4\\dM\\dD|5\\dM\\dD', giant.table.WT$cigar),])
    #gt.WT.del=cbind(gt.WT.del, oligos[match(gt.WT.del$oligo, rownames(oligos)),])

# merge in annotation information
    gt.WT.pe=cbind(gt.WT.p, oligos[match(gt.WT.p$oligo, oligos$oligo),])
    gt.NMD.pe=cbind(gt.NMD.p, oligos[match(gt.NMD.p$oligo, oligos$oligo),])

# more rearranging --------------------------------------------
    gWs=split(gt.WT.pe, gt.WT.pe$oligo)
    gNs=split(gt.NMD.pe, gt.NMD.pe$oligo)

    gtemp1=gt.WT.pe
    gtemp2=gt.NMD.pe
    colnames(gtemp2)=colnames(gt.WT.pe)
    gtemp1$expt='WT'
    gtemp2$expt='NMD'
    gBsA=rbind(gtemp1, gtemp2)
    gBs=split(gBsA, gBsA$oligo)
    gBsA=rbindlist(gBs) #do.call('rbind', gBs)
    rm(gtemp1)
    rm(gtemp2)
#-------------------------------------------------------------

# Fitting the oligo specific slope models  (note, these are slow)
#o.WT.fits=fit.oligo.model(gWs, t.var, colSums(gt.WT.pe[,9:13]),20)
#save(o.WT.fits, file = paste0(out.base.dir, 'processed/RData/o.WT.fits.RData2'))
load(paste0(out.base.dir, 'processed/RData/o.WT.fits.RData2'))
#o.N.fits=fit.oligo.model(gNs, t.var, colSums(gt.NMD.pe[,9:13]), 20)
#save(o.N.fits, file = paste0(out.base.dir, 'processed/RData/o.N.fits.RData2'))
load(paste0(out.base.dir, 'processed/RData/o.N.fits.RData2'))
#o.A.fits=fit.oligo.model.combined(gBs, t.var, colSums(gBsA[,9:13]),20)
#save(o.A.fits, file = paste0(out.base.dir, 'processed/RData/o.A.fits.RData2'))
load(paste0(out.base.dir, 'processed/RData/o.A.fits.RData2'))
#stopImplicitCluster()

# fit the big mixed model
o.cnt=sapply(o.A.fits, function(x) x$bc.cnt)
o.factor=rep(names(o.cnt), o.cnt)
o.gene=oligos$GENEID[match(names(o.cnt), rownames(oligos))]

o.cnt=rep((o.cnt), o.cnt)
o.intercept=sapply(o.A.fits, function(x) x$fixed.coeffs )
o.fit=do.call('rbind', o.intercept)
#merge it all together
raw.oligo.data=mapply(function(x,y) {
              x2= x[,c(9:13,14, 47, 120)]
              if(is.null(dim(y))) { y2=matrix(y,1,2) } else {y2=y }
              xx=cbind(x2, y2)
              binarized=ifelse(xx[,ncol(xx)]<(-.025), 0, 1)
              xx=cbind(xx,binarized)
                   }, x=gBs, y=o.intercept, SIMPLIFY=FALSE)
rawod.big=rbindlist(raw.oligo.data)

# Build and augment big.mm---------------------------------------------------------------
    big.mm=data.frame(gBsA, cnt=o.cnt, slope=o.fit[,2], p_intercept=o.fit[,1])
    big.mm$slope.binarized=ifelse(big.mm$slope>(-.025), 1,0)
    # For HMM
    big.mm$DA=as.factor(ifelse(big.mm$slope.binarized, '1', '2'))
    # now merge big.mm and oligo.stats ???
    #big.mm2=data.frame(oligo.stats[match(big.mm$oligo, oligo.stats$oligo),], big.mm, stringsAsFactors=F)
    #big.mm$close2end150=big.mm$dist_from_CDS_end
    #big.mm$close2end150[big.mm$close2end150>150]=NA

    #low complexity domain indicator
    big.mm$lcd=as.numeric(ifelse(big.mm$low.complexity.downstream.cnt>0, 1, 0))
    # categroize as having off target gRNAs
    big.mm$binot=ifelse(big.mm$cnt.offtarget==1,0,1)
    save(big.mm, file = paste0(out.base.dir, 'processed/RData/big.mm.RData'))
    #load(paste0(out.base.dir, 'processed/RData/big.mm.RData'))
#--------------------------------------------------------------------------------- 

rm(gBsA, gBs, gWs,gNs)

#merge information from oligo models
dW=t(sapply(o.WT.fits, function(x) { c(x$bc.cnt, x$feffs, x$p.value, unlist(x$cInt)  ) }))
dN=t(sapply(o.N.fits, function(x)  { c(x$bc.cnt,  x$feffs, x$p.value,  unlist(x$cInt)) }))
dA=t(sapply(o.A.fits, function(x)  { c(x$bc.cnt,  x$feffs[c(1,2)], x$p.value.time,  unlist(x$cInt), x$feffs[3], x$p.value.expt) }))
colnames(dW)=paste('WT', c('cnt', 'Intercept', 'slope', 'p.value', 'slope.025', 'slope.975'), sep='.')
colnames(dN)=paste('NMD',c('cnt', 'Intercept', 'slope', 'p.value', 'slope.025', 'slope.975'), sep='.')
colnames(dA)=paste('ALL',c('cnt', 'Intercept', 'slope', 'p.value', 'slope.025', 'slope.975', 'expt.eff', 'p.value.expt'), sep='.')

# Example mixture model
#mS=(Mclust(big.mm$slope, G=2, modelNames=c('E',  'V')))
#plot(mS, what='classification')

slope.mm=lmer(slope~cnt+p_intercept+expt+binot+(1|GENEID)+(1|oligo), data=big.mm, REML=TRUE)
slope.ranef=ranef(slope.mm, cond=TRUE)

gb=cbind(slope.ranef$GENEID[,1], as.vector(attr(slope.ranef$GENEID, 'postVar')))
colnames(gb)=c('slope.gene.blup', 'slope.gene.blup.postVar')
rownames(gb)=rownames(slope.ranef$GENEID)

binarized.mm=glmer(slope.binarized~cnt+p_intercept+expt+binot+(1|GENEID)+(1|oligo), family=binomial(link="logit"), data=big.mm, verbose=T)
binarized.ranef=ranef(binarized.mm, cond=TRUE)

gr=cbind(binarized.ranef$GENEID[,1], as.vector(attr(binarized.ranef$GENEID, 'postVar')))
colnames(gr)=c( 'binarized.gene.blup', 'binarized.gene.blup.postVar')
rownames(gr)=rownames(binarized.ranef$GENEID)

dA=cbind(dA, slope.ranef$oligo[,1], as.vector(attr(slope.ranef$oligo, 'postVar')), binarized.ranef$oligo[,1],as.vector(attr(binarized.ranef$oligo, 'postVar')) )
colnames(dA)[9:12]=c('slope.oligo.blup', 'slope.oligo.blup.postVar',
                     'binarized.oligo.blup', 'binarized.oligo.blup.postVar' )

# rewrite without using merge ()
os= data.frame(oligos, dW[match(oligos$oligo, rownames(dW)),], stringsAsFactors=F)
os= data.frame(os, dN[match(oligos$oligo, rownames(dN)),], stringsAsFactors=F)
os= data.frame(os, dA[match(oligos$oligo, rownames(dA)),], stringsAsFactors=F)
os= data.frame(os, gb[match(oligos$GENEID, rownames(gb)),], stringsAsFactors=F)
os= data.frame(os, gr[match(oligos$GENEID, rownames(gr)),], stringsAsFactors=F)
oligo.stats=os
rm(os)

# rename as oligos and resort 
oligo.stats=oligo.stats[order(oligo.stats$unique.Index),]
save(oligo.stats, file=paste0(out.base.dir, 'processed/RData/oligo.stats.1.RData'))
load(paste0(out.base.dir, 'processed/RData/oligo.stats.1.RData'))


#  HMM ================================================================================================
    # do HMM 
    hmm.out=doHMM(big.mm, oligo.stats)
    #save(hmm.out, file=paste0(out.base.dir, 'processed/RData/hmm.out.022217.RData'))
    load(paste0(out.base.dir, 'processed/RData/hmm.out.022217.RData'))

    # integrate data from hmm 2 state into oligo table 
    so=sapply(hmm.out, function(x) x$states.out)

    dso=data.frame(GENEID=rep(names(so), sapply(so, length )), 
                   dist_from_CDS_end=as.vector(unlist(sapply(so, names))),
                   hmm2=as.vector(unlist(so)))
    dso$gp=paste(dso$GENEID, dso$dist_from_CDS_end, sep=':')

    oligo.stats$gp=paste(oligo.stats$GENEID, oligo.stats$dist_from_CDS_end, sep=':')
    oligo.stats$hmm2=dso$hmm2[match(oligo.stats$gp, dso$gp)]

    # integrate data from hmm 3 state into oligo table
    so3=sapply(hmm.out, function(x) x$states.out3)
    dso3=data.frame(GENEID=rep(names(so3), sapply(so3, length )), 
                   dist_from_CDS_end=as.vector(unlist(sapply(so3, names))),
                   hmm3=as.vector(unlist(so3)))
    dso3$gp=paste(dso3$GENEID, dso3$dist_from_CDS_end, sep=':')
    oligo.stats$hmm3=dso3$hmm3[match(oligo.stats$gp, dso3$gp)]
    #save(oligo.stats, file=paste0(out.base.dir, 'processed/RData/oligo.stats.2.RData'))
    load(paste0(out.base.dir, 'processed/RData/oligo.stats.2.RData'))

    #-------------------------------------------------------------------
    plot.dir='/media/jbloom/d1/coupled_CRISPR/plots/HMM_v7/'
    makeHMMplots(plot.dir, oligo.stats, big.mm, conservation, dgsplit)
    
# ====================================================================================================

# confident alive oligos for example 
    (oligo.stats[which(-log10(oligo.stats$ALL.p.value)>5& 
                           oligo.stats$hmm2==1 & 
                           oligo.stats$binarized.oligo.blup>0.5 &
                           oligo.stats$hmm3==1)
                                ,])
    plot( oligo.stats$binarized.oligo.blup,  -log10(oligo.stats$ALL.p.value))



# Recalculate BLUP oligo and gene model and remove the essential genes 
formula.input=as.formula(slope.binarized~scale(cnt)+scale(p_intercept)+as.factor(expt)+scale(guide.GCcontent)+scale(PAMvariantCNT)+binot+(1|GENEID)+(1|oligo))
oligo.stats=doGLMER.oligo.gene(formula.input, data.in=big.mm[!big.mm$drop,], oligo.stats, lab='EssentialAll')

formula.input=as.formula(slope.binarized~scale(cnt)+scale(p_intercept)+scale(guide.GCcontent)+scale(PAMvariantCNT)+binot+(1|GENEID)+(1|oligo))
oligo.stats=doGLMER.oligo.gene(formula.input, data.in=big.mm[!big.mm$drop & (big.mm$expt == 'WT'),], oligo.stats, lab='EssentialWT')

formula.input=as.formula(slope.binarized~scale(cnt)+scale(p_intercept)+scale(guide.GCcontent)+scale(PAMvariantCNT)+binot+(1|GENEID)+(1|oligo))
oligo.stats=doGLMER.oligo.gene(formula.input, data.in=big.mm[!big.mm$drop & (big.mm$expt == 'NMD'),], oligo.stats, lab='EssentialNMD')
#============================================================================================
#save(oligo.stats, file=paste0(out.base.dir, 'processed/RData/oligo.stats.EssentialBlups.RData'))
load(paste0(out.base.dir, 'processed/RData/oligo.stats.EssentialBlups.RData'))

oligo.stats[oligo.stats$GENEID=='YIL051C',]
with(oligo.stats, {
par(mfrow=c(4,1))
hist(binarized.gene.blup, breaks=1000)
hist(binarized.gene.blup.EssentialAll,breaks=1000)
hist(binarized.gene.blup.EssentialWT,breaks=1000)
hist(binarized.gene.blup.EssentialNMD,breaks=1000)
                                })

# Generalized additive models to investigate effect of distance from end 
# Models  to get smoother lines SLOW
gamm=list()

gammALL     =gamm4(slope.binarized~s(dist_from_CDS_end) +scale(cnt) +as.factor(expt)+scale(PAMvariantCNT)+ scale(p_intercept)+ scale(guide.GCcontent)+binot, random=~(1|oligo), family=binomial(link='logit'), data=big.mm[!big.mm$drop,])
gammALL.dd  =gamm4(slope.binarized~s(dist_from_CDS_end) +scale(cnt) +as.factor(expt)+scale(PAMvariantCNT)+ scale(p_intercept)+ scale(guide.GCcontent)+binot, random=~(1|oligo), family=binomial(link='logit'), data=big.mm[!big.mm$drop & big.mm$domain.downstream,]) 
gammALL.ndd =gamm4(slope.binarized~s(dist_from_CDS_end) +scale(cnt) +as.factor(expt)+scale(PAMvariantCNT)+ scale(p_intercept)+ scale(guide.GCcontent)+binot, random=~(1|oligo), family=binomial(link='logit'), data=big.mm[!big.mm$drop & !big.mm$domain.downstream,])

gammWT      =gamm4(slope.binarized~s(dist_from_CDS_end) +scale(cnt)+ scale(p_intercept)+scale(PAMvariantCNT)+ scale(guide.GCcontent)+binot, random=~(1|oligo),family=binomial(link='logit'), data=big.mm[!big.mm$drop & (big.mm$expt == 'WT'),]) 
gammWT.dd   =gamm4(slope.binarized~s(dist_from_CDS_end) +scale(cnt)+ scale(p_intercept)+scale(PAMvariantCNT)+ scale(guide.GCcontent)+binot, random=~(1|oligo),family=binomial(link='logit'), data=big.mm[!big.mm$drop & (big.mm$expt == 'WT') & big.mm$domain.downstream,] )
gammWT.ndd  =gamm4(slope.binarized~s(dist_from_CDS_end) +scale(cnt)+ scale(p_intercept)+scale(PAMvariantCNT)+ scale(guide.GCcontent)+binot, random=~(1|oligo),family=binomial(link='logit'), data=big.mm[!big.mm$drop & (big.mm$expt == 'WT') & !big.mm$domain.downstream,] )

b3=big.mm[!big.mm$drop, ] # & (big.mm$expt == 'WT') , ]#& big.mm$dist_from_CDS_end<300,]
b3$domain.downstream=as.numeric(b3$domain.downstream)+1
b3$cnt=as.vector(scale(b3$cnt))
b3$p_intercept=as.vector(scale(b3$p_intercept))
b3$PAMvariantCNT=as.vector(scale(b3$PAMvariantCNT))
b3$guide.GCcontent=as.vector(scale(b3$guide.GCcontent))
b3$oligo=as.factor(b3$oligo)
#gammALL.split  =gamm4(slope.binarized~s(dist_from_CDS_end,by=domain.downstream)+cnt+p_intercept+PAMvariantCNT+guide.GCcontent, random=~(1|oligo),family=binomial(link='logit'), data=b3) 
#mcl=makeCluster()

gammALL.split  =bam(slope.binarized~s(dist_from_CDS_end,by=domain.downstream)+s(oligo, bs='re')+cnt+p_intercept+PAMvariantCNT+guide.GCcontent,family=binomial(link='logit'), data=b3, nthreads=50) 

df=data.frame(dist_from_CDS_end=c(1:300),
               domain.downstream=0,
               cnt=mean(b3$cnt),
               expt=0.5,
               PAMvariantCNT=mean(b3$PAMvariantCNT),
               p_intercept=mean(b3$p_intercept),
               guide.GCcontent=mean(b3$guide.GCcontent),
               oligo=1)

p1=predict(gammWT.split$gam, newdata=df, se=T)
df$domain.downstream=1
p2=predict(gammWT.split$gam, newdata=df, se=T)
plot(p1[[1]],ylim=c(-2,1), col='blue')
points(p2[[1]])

gammNMD     =gamm4(slope.binarized~s(dist_from_CDS_end) +scale(cnt)+ scale(p_intercept)+scale(PAMvariantCNT)+ scale(guide.GCcontent)+binot, random=~(1|oligo),family=binomial(link='logit'),  data=big.mm[!big.mm$drop & (big.mm$expt == 'NMD') ,]) 
gammNMD.dd  =gamm4(slope.binarized~s(dist_from_CDS_end) +scale(cnt)+ scale(p_intercept)+scale(PAMvariantCNT)+ scale(guide.GCcontent)+binot, random=~(1|oligo),family=binomial(link='logit'), data=big.mm[!big.mm$drop & (big.mm$expt == 'NMD') & big.mm$domain.downstream,])
gammNMD.ndd =gamm4(slope.binarized~s(dist_from_CDS_end) +scale(cnt)+ scale(p_intercept)+scale(PAMvariantCNT)+ scale(guide.GCcontent)+binot, random=~(1|oligo),family=binomial(link='logit'), data=big.mm[!big.mm$drop & (big.mm$expt == 'NMD') & !big.mm$domain.downstream,])
gammList=list(gammALL     =gammALL    , 
              gammWT      =gammWT     ,
              gammNMD     =gammNMD    ,
              gammALL.dd  =gammALL.dd ,
              gammALL.ndd =gammALL.ndd,
              gammWT.dd   =gammWT.dd  ,
              gammWT.ndd  =gammWT.ndd , 
              gammNMD.dd  =gammNMD.dd ,  
              gammNMD.ndd =gammNMD.ndd)
save(gammList, file= paste0(out.base.dir, 'processed/RData/gamm.fits.RData'))

#plot(g4$gam, xlim=rev(c(0,500)))
#points(big.mm3$dist_from_CDS_end, big.mm3$slope)
#optimizer="Nelder_Mead",
gALL=glmer(slope.binarized~ scale(cnt)+as.factor(expt)+scale(p_intercept)+scale(guide.GCcontent)+scale(dist_from_CDS_end)+scale(CDS_length)+
                            scale(PAMvariantCNT)+binot+lcd+evolvability+scale(mean.aa.perfect.conserved)+viable_annotation+scale(end.conservation)+as.numeric(domain.downstream)+
                            (1|GENEID)+(1|oligo), data=big.mm[!big.mm$drop &big.mm$dist_from_CDS_end<150,], family=binomial(link='logit'), control=glmerControl( optCtrl = list(maxfun = 1000000)),verbose=T)
#gALL2=allFit(gALL)
gALL.anova=(Anova(gALL, type='III'))
length(predict(gALL))
gALL.anova[order(gALL.anova[,1], decreasing=T),]

gWT=glmer(slope.binarized~ scale(cnt)+scale(p_intercept)+scale(guide.GCcontent)+scale(dist_from_CDS_end)+scale(CDS_length)+
                            scale(PAMvariantCNT)+binot+lcd+evolvability+scale(mean.aa.perfect.conserved)+ viable_annotation+scale(end.conservation)+as.factor(domain.downstream)+
                            (1|GENEID)+(1|oligo), data=big.mm[!big.mm$drop & (big.mm$dist_from_CDS_end<150) & big.mm$expt=='WT',], family=binomial(link='logit'), control=glmerControl(optCtrl = list(maxfun = 1000000)),verbose=T)
#gWT2=allFit(gWT)
gWT.anova=(Anova(gWT, type='III'))
length(predict(gWT))
gWT.anova[order(gWT.anova[,1], decreasing=T),]

gNMD=glmer(slope.binarized~ scale(cnt)+scale(p_intercept)+scale(guide.GCcontent)+scale(dist_from_CDS_end)+scale(CDS_length)+
                            scale(PAMvariantCNT)+binot+lcd+evolvability+scale(mean.aa.perfect.conserved)+ viable_annotation+scale(end.conservation)+as.factor(domain.downstream)+
                            (1|GENEID)+(1|oligo), data=big.mm[!big.mm$drop & (big.mm$dist_from_CDS_end<150) & big.mm$expt=='NMD',], family=binomial(link='logit'), control=glmerControl(optCtrl = list(maxfun = 1000000)),verbose=T)
#gNMD2=allFit(gNMD)
gNMD.anova=(Anova(gNMD, type='III'))
length(predict(gNMD))
gNMD.anova[order(gNMD.anova[,1], decreasing=T),]




relgrad <- with(gWT@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
####################################################################
# Figure 1 pilot testing 
git.dir='/media/jbloom/d1/coupled_CRISPR/Experiments/Preliminary/Pilot_092015/'
load(paste0(git.dir, 'summaryStats.RData'))
#load(paste0(git.dir, 'sampleTable.RData'))

summary.stats.table=t(do.call('rbind', summary.stats))
pilot.norm=t(t(summary.stats.table)/colSums(summary.stats.table))

#1 = matches haplotype no mismatches or indels in 50bp window of pam
#2 = indel in 50 bp window around pam
#3 = everything else
#4 = WT
par(oma=c(2,25,2,2))
pilot.drop=c(7,9:16, 17:19) #9:16, 23:25)
pilot.norm.drop=pilot.norm[,-pilot.drop]
pilot.norm.drop=pilot.norm.drop[,ncol(pilot.norm.drop):1]

bp=barplot(pilot.norm.drop, density=c(100,55,25,0),
           horiz=TRUE, xlim=c(0,1), yaxt='n', 
           beside=T, xlab='fraction of reads', space=c(0,3), width=2)
axis(2, bp[2,], colnames(pilot.norm.drop), las=2)
legend('bottomright',density=c(100,55,25,0) ,                      c('matches expected edit', 
                                                                        'indel near edit',
                                                                        'mismatch near edit', 
                                                                        'unedited'))


# Figure 2 slopes (Dubious vs essential)
par(oma=c(1,1,1,1))
with(big.mm, {
    plot(jitter(as.numeric(dubious)), slope, col='#00000010', xaxt='n', xlab='', ylab='barcode slope')
    axis(1,at=c(0,1), lab=c('dubious', 'essential'))  
    abline(h=-0.025, col='red')   
    text(c(.5,.5),c((-0.025-.02),(-0.025+.02)), c('dead', 'alive')) })

#-------------------------------------------------------------------------

# Figure 3 distance from end  (WT vs NMD) ---------------------------------------------------
gAp=plot(gammALL$gam)
gApD=plot(gammALL.dd$gam)
gApnD=plot(gammALL.ndd$gam)

dev.off()
plot(gAp[[1]]$x, gAp[[1]]$fit, type='l', xlim=c(200,0))
points(gAp[[1]]$x, gAp[[1]]$fit+gAp[[1]]$se, type='l', lty=2)
points(gAp[[1]]$x, gAp[[1]]$fit-gAp[[1]]$se, type='l', lty=2)


 plot_by_pos=with(oligo.stats[!oligo.stats$drop & !is.na(oligo.stats$binarized.oligo.blup.EssentialAll) & oligo.stats$dist_from_CDS_end<150,], {
   
    x=oligo.stats[!oligo.stats$drop & !is.na(oligo.stats$binarized.oligo.blup.EssentialAll) & oligo.stats$dist_from_CDS_end<150,]             
    blist=split(x$binarized.oligo.blup.EssentialAll[x$domain.downstream],x$dist_from_CDS_end[x$domain.downstream])
    blist2=split(x$binarized.oligo.blup.EssentialAll[!x$domain.downstream],x$dist_from_CDS_end[!x$domain.downstream])

    atx=as.numeric(names(blist))-.25
    atx2=as.numeric(names(blist2))+.25
    boxplot(blist,col='#ff000022',  border='#ff000022', xaxt='n',  boxwex=.25, at=atx, xlim=c(150,0))
    boxplot(blist2,col='#00000022',  border='#00000022', xaxt='n',  boxwex=.25, at=atx2, add=T)
    
    points(gApD[[1]]$x, gApD[[1]]$fit, type='l', lwd=2, col='red')
    points(gApD[[1]]$x, gApD[[1]]$fit+1.96*gApD[[1]]$se, type='l', lwd=1, lty=2, col='red')
    points(gApD[[1]]$x, gApD[[1]]$fit-1.96*gApD[[1]]$se, type='l', lwd=1, lty=2, col='red')
     
    points(gApnD[[1]]$x, gApnD[[1]]$fit, type='l', lwd=2, col='black')
    points(gApnD[[1]]$x, gApnD[[1]]$fit+1.96*gApnD[[1]]$se, type='l', lwd=1, lty=2, col='black')
    points(gApnD[[1]]$x, gApnD[[1]]$fit-1.96*gApnD[[1]]$se, type='l', lwd=1, lty=2, col='black')

    boxplot(binarized.oligo.blup.EssentialAll[domain.downstream]~dist_from_CDS_end[domain.downstream], col='#ff000022',  border='#ff000022', xaxs=F, xlim=c(150, 0), boxwex=.5, at=rev(atx))
    
    stripchart(binarized.oligo.blup.EssentialAll[domain.downstream]~dist_from_CDS_end[domain.downstream], vertical=T, add=T, pch=19, cex=.5, col='#ff000033')
    stripchart(binarized.oligo.blup.EssentialAll[!domain.downstream]~dist_from_CDS_end[!domain.downstream], vertical=T, add=T, pch=19, cex=.5, col='#00000033')
   
    #stripchart(binarized.oligo.blup.EssentialAll~dist_from_CDS_end, vertical=T, add=T, pch=19, cex=.5, col='#00000033')
    #stripchart(binarized.oligo.blup.EssentialAll~dist_from_CDS_end, vertical=T, add=T, pch=19, cex=.5, col='#00000033')

    # stripchart(     binarized.oligo.blup~dist_from_CDS_end, vertical=T, add=T, pch=19, cex=.5)
    #boxplot(binarized.oligo.blup[domain.downstream]~dist_from_CDS_end[domain.downstream], col='#ff000033',  border='#ff000033', xaxs=F, xlim=c(150, 0))
    #boxplot(binarized.oligo.blup[!domain.downstream]~dist_from_CDS_end[!domain.downstream], col='#00000033',  border='#00000033', xaxs=F, xlim=c(150, 0),  xlim=range(dist_from_CDS_end[domain.downstream]) ,add=T )                                                                                                                                                      
                                                                                                                                                                      
     #points(smooth.spline(dist_from_CDS_end, binarized.oligo.blup), type='l', col='blue', lwd=2)
     #points(smooth.spline(dist_from_CDS_end[domain.downstream], binarized.oligo.blup[domain.downstream]), type='l', col='red', lwd=2)
     #points(smooth.spline(dist_from_CDS_end[!domain.downstream], binarized.oligo.blup[!domain.downstream]), type='l', col='black', lwd=2)
           }
)

# Figure 4 domains / other features of models


# Figure 5 gene specific analysis  (GO)
idx=match(unique(oligo.stats$GENEID), oligo.stats$GENEID) 
stripchart(oligo.stats$binarized.gene.blup[idx]~oligo.stats$dubious[idx], vertical=T, xlab='DUBIOUS ORF', main='Gene Blups',
           ylab='gene persisitence score (+ = persisting)', jitter=.1, method='jitter', pch=21)

            
# WT vs NMD
# WT vs NMD as Figure 3 (oligo blups by position)
##
pv <- predict(gammALL$gam,big.mm[!big.mm$drop,], type="response",se=TRUE) 
plot(big.mm[!big.mm$drop,]$dist_from_CDS_end, pv$fit)

lines(mydata$x0,pv$fit) 
lines(mydata$x0,pv$fit+2*pv$se.fit,lty=2) 
lines(mydata$x0,pv$fit-2*pv$se.fit,lty=2) 





# scale(H0_w)+
#scale(CDS_length)+
#binot+
# subset of data out to 125
# 2 step ... fit full model... then remove terms that explain nothin 
slopeA=glmer(slope.binarized~
                scale(cnt)+
                as.factor(expt)+
                scale(p_intercept)+
                scale(guide.GCcontent)+
                scale(close2end150)+
                scale(CDS_length)+
                PAMvariantCNT+
                binot+
                lcd+
                evolvability+
                mean.aa.perfect.conserved+
                viable_annotation+
                end.conservation+
                as.numeric(domain.downstream)+-1+
                (1|GENEID)+(1|oligo), data=big.mm3, 
            family=binomial(link="logit"),
            control=glmerControl(optimizer="Nelder_Mead", optCtrl = list(maxfun = 1e6)) )
            
            REML=TRUE)
    #summary(slope5)
    at=(anova(slopeA, test='F'))
    length(predict(slopeA))
    at[order(at[,1], decreasing=T),]


    
    
    
    
    
    
    
    plot_by_pos=with(oligo.stats[!oligo.stats$drop & !is.na(oligo.stats$binarized.oligo.blup) & oligo.stats$dist_from_CDS_end<150,], {
    #boxplot((exp(binarized.oligo.blup)/(1+exp(binarized.oligo.blup)))~dist_from_CDS_end, col='#ff000033',  border='#ff000033', xaxs=F, xlim=c(100, 0), ylim=c(0,1))
    #stripchart((exp(binarized.oligo.blup)/(1+exp(binarized.oligo.blup))) ~dist_from_CDS_end, vertical=T, add=T, pch=19, cex=.5)
    #boxplot(exp(binarized.oligo.blupR.WT)~dist_from_CDS_end, col='#ff000033',  border='#ff000033', xaxs=F, xlim=c(100, 0), ylim=c(0,20))
    #stripchart(exp(binarized.oligo.blupR.WT) ~dist_from_CDS_end, vertical=T, add=T, pch=19, cex=.5)
   
    boxplot(binarized.oligo.blupR.WT~dist_from_CDS_end, col='#ff000033',  border='#ff000033', xaxs=F, xlim=c(150, 0))
    stripchart(binarized.oligo.blupR.WT~dist_from_CDS_end, vertical=T, add=T, pch=19, cex=.5)

    # stripchart(     binarized.oligo.blup~dist_from_CDS_end, vertical=T, add=T, pch=19, cex=.5)
    #boxplot(binarized.oligo.blup[domain.downstream]~dist_from_CDS_end[domain.downstream], col='#ff000033',  border='#ff000033', xaxs=F, xlim=c(150, 0))
    #boxplot(binarized.oligo.blup[!domain.downstream]~dist_from_CDS_end[!domain.downstream], col='#00000033',  border='#00000033', xaxs=F, xlim=c(150, 0),  xlim=range(dist_from_CDS_end[domain.downstream]) ,add=T )                                                                                                                                                      
                                                                                                                                                                      
     #points(smooth.spline(dist_from_CDS_end, binarized.oligo.blup), type='l', col='blue', lwd=2)
     #points(smooth.spline(dist_from_CDS_end[domain.downstream], binarized.oligo.blup[domain.downstream]), type='l', col='red', lwd=2)
     #points(smooth.spline(dist_from_CDS_end[!domain.downstream], binarized.oligo.blup[!domain.downstream]), type='l', col='black', lwd=2)
           }
)

# Potential Figure, essentiality by GENE including non-stops
idx=match(unique(oligo.stats$GENEID), oligo.stats$GENEID) 
stripchart(oligo.stats$binarized.gene.blup[idx]~oligo.stats$dubious[idx], vertical=T, xlab='DUBIOUS ORF', main='Gene Blups',
           ylab='gene persisitence score (+ = persisting)', jitter=.1, method='jitter', pch=21)

wilcox.test(oligo.stats$binarized.gene.blup[idx]~oligo.stats$dubious[idx])$p.value

# end of gene conservation specific signature ?

# GO Enrichment analysis 
# extract oligo gene blups
testGenes=oligo.stats$binarized.gene.blup[match(unique(oligo.stats[!oligo.stats$drop,]$GENEID), oligo.stats$GENEID)]
names(testGenes)=oligo.stats$GENEID[match(unique(oligo.stats[!oligo.stats$drop,]$GENEID), oligo.stats$GENEID)]
testGenes=sort(testGenes, decreasing=T)
GO.out=doGO(testGenes)
head(GO.out[[1]])

fdr.out=qvalue(unlist(sapply(GO.out, function(x) x$p.value)), lambda=seq(.3,.5,.05), fdr.level=.01)
max(fdr.out$pvalues[fdr.out$significant])
#0.00086
thisOntology='CC'
GOData = new("topGOdata", ontology=thisOntology, allGenes = testGenes, geneSel=geneSel.fx, annot = annFUN.gene2GO, gene2GO = gene2GOList, nodeSize=5)
goID='GO:0005681'
plot(showGroupDensity(GOData, goID, ranks=T))
title('spliceosomal complex' )

#genesinterms=genesInTerm(GOData, goID)
WriteXLS(GO.out, '/home/jbloom/Dropbox/Public/GO_gene_blups_022217.xls', SheetNames=names(GO.out))






#e^x / (1 + e^x)
plot(t.var, c(0,0,0,0,0), ylim=c(0,500))
with(big.mm, 
     for(i in 1:nrow(big.mm)){
     points(t.var, c(WT_0_rpkm[i],  WT_24_rpkm[i]  ,  WT_48_rpkm[i], WT_72_rpkm[i]  ,  WT_96_rpkm[i]), type='l', col='#00000011')
     })
X11()
plot(t.var, c(0,0,0,0,0), ylim=c(0,500))
with(big.mm[big.mm$dubious,], 
     for(i in 1:nrow(big.mm[big.mm$dubious,])){
     points(t.var, c(WT_0_rpkm[i],  WT_24_rpkm[i]  ,  WT_48_rpkm[i], WT_72_rpkm[i]  ,  WT_96_rpkm[i]), type='l', col='#00000011')
     })



#$  WT_0_rpkm                    : num  48.2 46.9 48.2 143.3 164 ...
#$ WT_24_rpkm                   : num  13.2 275 37.7 909.5 75.4 ...
#$ WT_48_rpkm                   : num  0 642 0 2452 0 ...
#$ WT_72_rpkm                   : num  0 375.82 0 1573.08 1.07 ...
#$ WT_96_rpkm                   : num  0 285 0 2342 0 ...


dff=na.omit(data.frame(oligo.stats$dist_from_CDS_end, oligo.stats$binarized.oligo.blupR))
dff[,2]=dff[,2]+min(dff[,2])
dff3=dff[dff[,1]<200,]

summary(nlsfit(dff3, model=6, start=c(a=-1,b=-1)))

with(oligo.stats[!oligo.stats$drop & !is.na(oligo.stats$binarized.oligo.blup) & oligo.stats$dist_from_CDS_end<150,], {
           
         summary(nlsfit(binarized.oligo.blup~dist_from_CDS_end, model=1))
})
dm=Mclust(dff[,2], G=2, modelNames=c('E',  'E'))

sapply(split(log10(dm$z[,2]/dm$z[,1]), dff[,1]), sum)
plot(sapply(split(log10(dm$z[,2]/dm$z[,1]), dff[,1]), mean) ,xlim=c(0,200), ylim=c(-2,2) )#, xlim=c(0,200))
plot(sapply(split(log10(dm$z[,2]/dm$z[,1]), dff[,1]), length) ,xlim=c(0,200) )#, xlim=c(0,200))

slopeAR=glmer(slope.binarized~
            scale(cnt)+
            as.factor(expt)+
            scale(p_intercept)+
            scale(guide.GCcontent)+
            (1|close2end50)+
            scale(PAMvariantCNT)+-1+
            (1|GENEID)+
            (1|oligo), data=big.mm3, family=binomial(link="logit"),control=glmerControl(optimizer="Nelder_Mead", optCtrl = list(maxfun = 500000)))#REML=TRUE)

getInitial(slope.binarized ~ SSlogis(dist_from_CDS_end, Asym, xmid,scal), data = big.mm3)

initVals <- getInitial(slope.binarized ~ SSasympOff(dist_from_CDS_end), data = b)

slopeAR=glmer(slope.binarized~
            scale(cnt)+
            as.factor(expt)+
            scale(p_intercept)+
            scale(guide.GCcontent)+
            SS+
            scale(PAMvariantCNT)+-1+
            (1|GENEID)+
            (1|oligo), data=big.mm3, family=binomial(link="logit"),control=glmerControl(optimizer="Nelder_Mead", optCtrl = list(maxfun = 500000)))#REML=TRUE)

b=big.mm3
b=b[!is.na(b$slope.binarized) & !is.na(b$close2end100),]


par(mfrow=c(2,1))
#plot(g4$gam, xlim=rev(c(0,500)), main='both')
plot(g5$gam, xlim=rev(c(0,500)), main='WT', ylim=c(-1,1))
abline(v=100)
 abline(h=0)

plot(g6$gam, xlim=rev(c(0,500)), main='NMD', ylim=c(-1,1))
abline(v=100)
 abline(h=0)

     
     # #xlim=rev(c(0,500)))




rg=ranef(g4$mer)
plot(oligo.stats$dist_from_CDS_end[match(rownames(rg$oligo), oligo.stats$oligo)], rg$oligo[,1], xlim=rev(c(0,500))) #b$dist_from_CDS_end,, xlim=rev(c(0,500)))
boxplot(rg$oligo[,1]~ oligo.stats$dist_from_CDS_end[match(rownames(rg$oligo), oligo.stats$oligo)] , xlim=rev(c(0,100)) )

, rg$oligo[,1], xlim=rev(c(0,500))) #b$dist_from_CDS_end,, xlim=rev(c(0,500)))



getInitial(slope.binarized ~ SSlogis(dist_from_CDS_end, Asym, xmid,scal), data = b)

pslope=predict(slopeAR)
ra=ranef(slopeAR)
boxplot(ra$oligo[,1]~oligo.stats$dist_from_CDS_end[match(rownames(ra$oligo), oligo.stats$oligo)])

rslopeAR=residuals(slopeAR)
big.mm3$rslope=rep(NA, nrow(big.mm3))
big.mm3$rslope[match(names(rslopeAR), rownames(big.mm3))]=as.vector(rslopeAR)
slopeAR2=lmer(rslope~
            scale(cnt)+
            as.factor(expt)+
            scale(p_intercept)+
            scale(guide.GCcontent)+
            scale(PAMvariantCNT)+
            +dist_from_CDS_end
            (1|GENEID)+
            (1|oligo)+ 
            evolvability+
            mean.aa.perfect.conserved+
            viable_annotation+
            end.conservation+
            as.numeric(domain.downstream)+-1,
        data=big.mm3) #, family=binomial(link="logit"),control=glmerControl(optimizer="Nelder_Mead", optCtrl = list(maxfun = 500000)))#REML=TRUE)
    at2=(anova(slopeAR2, test='F'))





#summary(slope5)
    at=(anova(slopeA, test='F'))
    length(predict(slopeA))
    at[order(at[,1], decreasing=T),]

 lcd+
            evolvability+
            mean.aa.perfect.conserved+
            viable_annotation+
            end.conservation+
            as.numeric(domain.downstream)+-1+



    # scale(H0_w)+
    #scale(CDS_length)+
    #binot+
    slopeA=lmer(slope~
            scale(cnt)+
            as.factor(expt)+
            scale(p_intercept)+
            scale(guide.GCcontent)+
            (1|close2end100)+
            PAMvariantCNT+
            lcd+
            evolvability+
            mean.aa.perfect.conserved+
            viable_annotation+
            end.conservation+
            as.numeric(domain.downstream)+-1+
             (1|GENEID)+(1|oligo), data=big.mm3, REML=TRUE)
    #summary(slope5)
    at=(anova(slopeA, test='F'))
    length(predict(slopeA))
    at[order(at[,1], decreasing=T),]






    slopeA.100=lmer(slope~
            scale(cnt)+
            as.factor(expt)+
            scale(p_intercept)+
            scale(guide.GCcontent)+
            scale(close2end100)+
            PAMvariantCNT+
            lcd+
            evolvability+
            mean.aa.perfect.conserved+
            viable_annotation+
            end.conservation+
            as.numeric(domain.downstream)+-1+
             (1|GENEID)+(1|oligo), data=big.mm3, REML=FALSE)
    #summary(slope5)
    at.100=(anova(slopeA.100, test='F'))
    length(predict(slopeA.100))
    at.100[order(at.100[,1], decreasing=T),]

    slopeA.100.int=lmer(slope~
            scale(cnt)+
            scale(p_intercept)+
            scale(guide.GCcontent)+
            scale(close2end100)*as.factor(expt)+
            PAMvariantCNT+
            lcd+
            evolvability+
            mean.aa.perfect.conserved+
            viable_annotation+
            end.conservation+
            as.numeric(domain.downstream)+-1+
             (1|GENEID)+(1|oligo), data=big.mm3, REML=FALSE)
    #summary(slope5)
    at.100.int=(anova(slopeA.100.int, test='F'))
    length(predict(slopeA.100.int))
    at.100.int[order(at.100.int[,1], decreasing=T),]
   

slopeA.100.intG=glmer(slope.binarized~
                    scale(cnt)+
            scale(p_intercept)+
            scale(guide.GCcontent)+
            scale(close2end100)*as.factor(expt)+
            PAMvariantCNT+
            lcd+
            evolvability+
            mean.aa.perfect.conserved+
            viable_annotation+
            end.conservation+
            as.numeric(domain.downstream)+-1+
             (1|GENEID)+(1|oligo), data=big.mm3, family=binomial(link="logit"),control=glmerControl(optimizer="Nelder_Mead")) # REML=FALSE)
    #summary(slope5)
    at.100.intG=(anova(slopeA.100.intG, test='F'))
    length(predict(slopeA.100.int))
    at.100.int[order(at.100.int[,1], decreasing=T),]







gc2=oligo.stats$guide.GCcontent
gc2[gc2<.25]=NA
gc2[gc2>.6]=NA
cor.test(gc2, oligo.stats$binarized.oligo.blup, method='pearson')

cor.test(oligo.stats$guide.GCcontent, oligo.stats$binarized.oligo.blup, method='spearman')

stripchart(oligo.stats$binarized.oligo.blup~oligo.stats$guide.GCcontent, vertical=T, method='jitter',col=ifelse(oligo.stats$dubious, 'red', 'black'))

boxplot(oligo.stats$binarized.oligo.blup~oligo.stats$guide.GCcontent, add=T)




plot(jitter(oligo.stats$guide.GCcontent), oligo.stats$binarized.oligo.blup, col=ifelse(oligo.stats$dubious, 'red', 'black'))
boxplot(oligo.stats$binarized.oligo.blup~oligo.stats$guide.GCcontent, add=T)



# additional info 82, 116
idx=match(unique(oligo.stats$GENEID), oligo.stats$GENEID) 
with(oligo.stats[ seq(1:nrow(oligo.stats)) %in% idx & !oligo.stats$drop & !is.na(oligo.stats$binarized.oligo.blup), ], {
         plot( log2(YEPD.mean), binarized.gene.blup,  ylab='gene persistence', xlab='log2(GFP)')
         cr1=cor.test( log2(YEPD.mean), binarized.gene.blup, method='spearman') #,  ylab='gene persistence', xlab='log2(GFP)')
        legend('topright', paste0('rho=', round(cr1$estimate,2), '      p< ', format.pval(cr1$p.value)))
                       })

# plot by conservation (gene level)
with(oligo.stats[ seq(1:nrow(oligo.stats)) %in% idx & !oligo.stats$drop & !is.na(oligo.stats$binarized.oligo.blup), ], {
    plot(y=binarized.gene.blup,x=mean.aa.perfect.conserved, ylab='gene persistence', xlab='conservation score', xlim=c(0.6,1))
    cr1=cor.test(mean.aa.perfect.conserved,binarized.gene.blup, method='spearman')
    legend('topleft', paste0('rho=', round(cr1$estimate,2), '      p< ', format.pval(cr1$p.value)))
                       })

# half life
#plot(y=red.effs$GENEID[,1], x=log2(argdf$Corrected.Half.Life), ylab='gene persistence', xlab='log2( protein half-life)')
#cr1=cor.test(log2(argdf$Corrected.Half.Life), red.effs$GENEID[,1], method='spearman')
#legend('topright', paste0('rho=', round(cr1$estimate,2), '      p< ', format.pval(cr1$p.value)))

# plot by dn_ds (H0_w) 

with(oligo.stats[ seq(1:nrow(oligo.stats)) %in% idx & !oligo.stats$drop & !is.na(oligo.stats$binarized.oligo.blup), ], {
    plot(y=binarized.gene.blup,x=H0_w, ylab='gene persistence', xlab='dN/dS', xlim=c(0,.35))
    cr1=cor.test(H0_w, binarized.gene.blup, method='spearman')
    legend('topleft', paste0('rho=', round(cr1$estimate,2), '      p< ', format.pval(cr1$p.value)))
                       })

# haploinsufficiency
with(oligo.stats[ seq(1:nrow(oligo.stats)) %in% idx & !oligo.stats$drop & !is.na(oligo.stats$binarized.oligo.blup), ], {
    plot(y=binarized.gene.blup,x=HET_AV, ylab='gene persistence', xlab='haploinsufficiency')
    cr1=cor.test(HET_AV, binarized.gene.blup, method='spearman')
    legend('topleft', paste0('rho=', round(cr1$estimate,2), '      p< ', format.pval(cr1$p.value)))
    cr2=cor.test(as.numeric(as.factor(slow_ypd_het)), binarized.gene.blup, method='spearman')
    print(cr2)
                       })
# length
with(oligo.stats[ seq(1:nrow(oligo.stats)) %in% idx & !oligo.stats$drop & !is.na(oligo.stats$binarized.oligo.blup), ], {
    plot(y=binarized.gene.blup,x=CDS_length, ylab='gene persistence', xlab='CDS length')
    cr1=cor.test(CDS_length, binarized.gene.blup, method='spearman')
    legend('topleft', paste0('rho=', round(cr1$estimate,2), '      p< ', format.pval(cr1$p.value)))
                       })

# yeast annotated viability info
# not restriciting to be only annotated in s288c or restricting on type of data

with(oligo.stats[ seq(1:nrow(oligo.stats)) %in% idx & !oligo.stats$drop & !is.na(oligo.stats$binarized.oligo.blup), ], {
stripchart(binarized.gene.blup~viable_annotation, vertical=TRUE, pch=21,method='jitter', xlim=c(0.5,2.5), ylab='gene persistence', xlab='any SGD annotation for viable',
           sub=paste0('wilcox p< ', format.pval(wilcox.test(binarized.gene.blup~viable_annotation)$p.value )))
                       
        boxplot(binarized.gene.blup~viable_annotation,add=T,xaxt='n')
                            }
)

# evolvability 
with(oligo.stats[ seq(1:nrow(oligo.stats)) %in% idx & !oligo.stats$drop & !is.na(oligo.stats$binarized.oligo.blup), ], {
stripchart(binarized.gene.blup~evolvability, vertical=T, method='jitter', pch=21, xlim=c(0.5,4.5), ylab='gene persistence')
boxplot(binarized.gene.blup~evolvability, add=T, xaxt='n') 
pairwise.t.test(binarized.gene.blup,evolvability, p.adjust='bonferroni')      
})


# human complements null
with(oligo.stats[ seq(1:nrow(oligo.stats)) %in% idx & !oligo.stats$drop & !is.na(oligo.stats$binarized.oligo.blup), ], {
    stripchart(binarized.gene.blup~human_complements_null-1, vertical=T, method='jitter', pch=21, xlim=c(0.5,2.5), ylab='gene persistence', xlab='human CDS complements null')
    boxplot(binarized.gene.blup~human_complements_null, add=T, xaxt='n')
    wilcox.test(binarized.gene.blup~human_complements_null, add=T, xaxt='n')
})







dim(big.mm3)

R> sum(is.na(big.mm3$mean.aa.perfect.conserved))
[1] 8335
R> sum(is.na(big.mm3$evolvability))
[1] 3453
R> sum(is.na(big.mm3$HET_AV))
[1] 4470
R> sum(is.na(big.mm3$viable_annotation))
[1] 0
R> sum(is.na(big.mm3$evolvability))
[1] 3453
R> sum(is.na(big.mm3$human_complements_null))
[1] 49978
R> sum(is.na(big.mm3$YEPD.mean))
[1] 37708
R> sum(is.na(big.mm3$end.conservation))
[1] 8463



slope4=lmer(slope~scale(cnt)+
            as.factor(expt)+
            scale(p_intercept)+
            scale(dist_from_CDS_end)+
            scale(H0_w)+
            scale(CDS_length)+
            +lcd+evolvability+mean.aa.perfect.conserved+viable_annotation+human_complements_null+scale(end.conservation)+scale(HET_AV)+
            as.numeric(domain.downstream)+scale(log2(YEPD.mean))-1+
             (1|GENEID)+(1|oligo), data=big.mm3, REML=FALSE)
#summary(slope4)
(anova(slope4, test='F'))
length(predict(slope4))
15192

# Main model here
binot=ifelse(big.mm3$cnt.offtarget==1,0,1)

slope5=lmer(slope~scale(cnt)+
            as.factor(expt)+
            scale(p_intercept)+
            scale(dist_from_CDS_end)+
            scale(H0_w)+
            scale(CDS_length)+
            scale(guide.GCcontent)+
            binot+
            PAMvariantCNT+
            +lcd+evolvability+mean.aa.perfect.conserved+viable_annotation+end.conservation+
            as.numeric(domain.downstream)+-1+
             (1|GENEID)+(1|oligo), data=big.mm3, REML=FALSE)
#summary(slope5)
at=(anova(slope5, test='F'))
length(predict(slope5))
#] 68087




# not close to end 
slopeN=lmer(slope~scale(cnt)+
            as.factor(expt)+
            scale(p_intercept)+
            scale(notclose2end)+
            scale(H0_w)+
            scale(CDS_length)+
            +lcd+evolvability+mean.aa.perfect.conserved+viable_annotation+end.conservation+
            as.numeric(domain.downstream)+-1+(1|GENEID)+(1|oligo), data=big.mm3, REML=FALSE)
#summary(slope5)
(anova(slopeN, test='F'))
length(predict(slopeN))

# interaction
slope16=lmer(slope~scale(cnt)+
            (expt)*scale(close2end100)+
            scale(p_intercept)+
            scale(H0_w)+
            scale(CDS_length)+
            +lcd+evolvability+mean.aa.perfect.conserved+viable_annotation+end.conservation+
            as.numeric(domain.downstream)+-1+
             (1|GENEID)+(1|oligo), data=big.mm3, REML=FALSE)
(anova(slope16, test='F'))

big.mm4=big.mm3
big.mm4=big.mm4[big.mm4$expt=='WT',]

big.mm5=big.mm3
big.mm5=big.mm5[big.mm5$expt=='NMD',]

# WT only
slopeWT=lmer(slope~scale(cnt)+
            scale(close2end100)+
            scale(p_intercept)+
            scale(H0_w)+
            scale(CDS_length)+
            +lcd+evolvability+mean.aa.perfect.conserved+viable_annotation+end.conservation+
            as.numeric(domain.downstream)+-1+
             (1|GENEID)+(1|oligo), data=big.mm4, REML=FALSE)
(anova(slopeWT, test='F'))


slopeNMD=lmer(slope~scale(cnt)+
            scale(close2end100)+
            scale(p_intercept)+
            scale(H0_w)+
            scale(CDS_length)+
            +lcd+evolvability+mean.aa.perfect.conserved+viable_annotation+end.conservation+
            as.numeric(domain.downstream)+-1+
             (1|GENEID)+(1|oligo), data=big.mm5, REML=FALSE)

(anova(slopeNMD, test='F'))



slope6=lmer(slope~scale(cnt)+
            as.factor(expt)+
            scale(p_intercept)+
            scale(close2end100)+
            scale(H0_w)+
            scale(CDS_length)+
            +lcd+evolvability+mean.aa.perfect.conserved+viable_annotation+end.conservation+
            as.numeric(domain.downstream)+-1+
             (1|GENEID)+(1|oligo), data=big.mm3, REML=FALSE)
#summary(slope6)
(anova(slope6, test='F'))
length(predict(slope6))

#[1] 30234
slope7=lmer(slope~scale(cnt)+
            as.factor(expt)+
            scale(p_intercept)+
            scale(close2end50)+
            scale(H0_w)+
            scale(CDS_length)+
            +lcd+evolvability+mean.aa.perfect.conserved+viable_annotation+end.conservation+
            as.numeric(domain.downstream)+-1+
             (1|GENEID)+(1|oligo), data=big.mm3, REML=FALSE)
#summary(slope5)
(anova(slope7, test='F'))
length(predict(slope7))



slope8=lmer(slope~scale(cnt)+
            scale(p_intercept)+
            scale(close2end50)+
            scale(H0_w)+
            scale(CDS_length)+
            +lcd
            +mean.aa.perfect.conserved
            +end.conservation+
            as.numeric(domain.downstream)+-1+
             (1|GENEID)+(1|oligo)+(1|viable_annotation)+(1|evolvability)+(1|expt), data=big.mm3, REML=FALSE)
(anova(slope8, test='F'))












lcd=as.numeric(ifelse(big.mm2$low.complexity.downstream.cnt>0, 1, 0))
slope2=lmer(slope~scale(cnt)+scale(p_intercept)+scale(close2end)+scale(H0_w)+scale(CDS_length)+lcd+
            as.numeric(domain.downstream)+scale(end.conservation)-1+
            scale(log2(YEPD.mean))+(1|GENEID)+(1|oligo)+(1|expt), data=big.mm2, REML=FALSE, na.action=na.pass)
summary(slope2)
(anova(slope2, test='F'))
 test=predict(slope2)



slopeg=glmer(slope.binarized~scale(cnt)+scale(p_intercept)+scale(close2end)+scale(H0_w)+scale(CDS_length)+lcd+
            as.numeric(domain.downstream)+scale(end.conservation)-1+
            scale(log2(YEPD.mean))+(1|GENEID)+(1|oligo), data=big.mm2,family=binomial(link="logit"),glmerControl(optimizer="bobyqa",optCtrl = list(maxfun = 500000)) )
summary(slopeg)

slope2=lmer(slope~scale(cnt)+scale(p_intercept)+scale(close2end)+H0_w+CDS_length+mean.aa.perfect.conserved+-1+log2(YEPD.mean)+(1|GENEID)+(1|oligo)+(1|expt), data=big.mm2, REML=FALSE)
summary(slope2)
(anova(slope2, test='F'))


b2=glmer(slope.binarized~cnt+p_intercept+expt+H0_w+HET_AV+CDS_length+mean.aa.perfect.conserved+evolvability+(1|GENEID)+(1|oligo),
         family=binomial(link="logit"), data=big.mm)

(anova(b2, test='F'))

test=sapply(split(lcd, big.mm2$dist_from_CDS_end), function(x) sum(x)/length(x))


slope.mm=lmer(slope~cnt+p_intercept+expt+(1|GENEID)+(1|oligo), data=big.mm, REML=TRUE)
slope.ranef=ranef(slope.mm, cond=TRUE)

gb=cbind(slope.ranef$GENEID[,1], as.vector(attr(slope.ranef$GENEID, 'postVar')))
colnames(gb)=c('slope.gene.blup', 'slope.gene.blup.postVar')
rownames(gb)=rownames(slope.ranef$GENEID)

binarized.mm=glmer(slope.binarized~cnt+p_intercept+expt+(1|GENEID)+(1|oligo), family=binomial(link="logit"), data=big.mm)
binarized.ranef=ranef(binarized.mm, cond=TRUE)


os= data.frame(os, dN[match(oligos$oligo, rownames(dN)),], stringsAsFactors=F)
















