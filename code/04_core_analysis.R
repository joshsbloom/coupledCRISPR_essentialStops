library(rbamtools)
library(parallel)
library(lme4)
library(stringdist)
library(mgcv)
library(bio3d)
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
library(optimx)
library(car)
#library(msm)
#library(afex)

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


source('/media/jbloom/d1/coupled_CRISPR/code/mobsHMM.R')

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
    par(mfrow=c(1,1))
    barplot(rbind(expt.counts.red2, expt.counts.red-expt.counts.red2, expt.counts-expt.counts.red),
            ylab='total reads', xlab='sample', main='reads per sample', col=c('black', 'grey', 'white'))
    legend('topright', c( paste0(round(sum(expt.counts)/1e6,1), ' x 1e6 total reads'),
                          paste0(round(sum(expt.counts.red)/1e6,1), ' x 1e6 total reads\n (reads per barcode filter)')  ,
                          paste0(round(sum(expt.counts.red2)/1e6,1), ' x 1e6 total reads\n (reads per barcode and perfect match repair)')),
                fill=c('black', 'grey', 'white'),
                   cex=1) 
   dev.off() 
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
    load(paste0(out.base.dir, 'processed/RData/big.mm.RData'))
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

os= data.frame(oligos, dW[match(oligos$oligo, rownames(dW)),], stringsAsFactors=F)
os= data.frame(os, dN[match(oligos$oligo, rownames(dN)),], stringsAsFactors=F)
os= data.frame(os, dA[match(oligos$oligo, rownames(dA)),], stringsAsFactors=F)

#with slopes and gene
formula.input=as.formula(slope~cnt+p_intercept+expt+binot+U6.Terminator+Score+guide.GCcontent+PAMvariantCNT+(1|GENEID)+(1|oligo))
oligo.stats=doLMER.oligo.gene(formula.input, data.in=big.mm, os, lab='')
#with slopes, no gene
formula.input=as.formula(slope~cnt+p_intercept+expt+binot+U6.Terminator+Score+guide.GCcontent+PAMvariantCNT+(1|oligo))
oligo.stats=doLMER.oligo.gene(formula.input, data.in=big.mm, oligo.stats, lab='ng', doGene=FALSE)

formula.input=as.formula(slope.binarized~cnt+p_intercept+expt+binot+U6.Terminator+Score+guide.GCcontent+PAMvariantCNT+(1|GENEID)+(1|oligo))
oligo.stats=doGLMER.oligo.gene(formula.input, data.in=big.mm, oligo.stats, lab='')
#with slopes, no gene
formula.input=as.formula(slope.binarized~cnt+p_intercept+expt+binot+U6.Terminator+Score+guide.GCcontent+PAMvariantCNT+(1|oligo))
oligo.stats=doGLMER.oligo.gene(formula.input, data.in=big.mm, oligo.stats, lab='ng', doGene=FALSE)
oligo.stats=oligo.stats[order(oligo.stats$unique.Index),]

save(oligo.stats, file=paste0(out.base.dir, 'processed/RData/oligo.stats.1.RData'))
load(paste0(out.base.dir, 'processed/RData/oligo.stats.1.RData'))

#save.image(paste0(out.base.dir, 'processed/RData/workspace03092017.RData'))

#  HMM ================================================================================================
    # do HMM 
    hmm.out=doHMM(big.mm, oligo.stats)
    oligo.stats = data.frame(oligo.stats, hmm.out[match(oligos$oligo, hmm.out$name),])
 
    acnt=sapply(split(oligo.stats$hmm2, oligo.stats$GENEID), function(x) sum(x=='A',na.rm=T))
    afrac=sapply(split(oligo.stats$hmm2, oligo.stats$GENEID), function(x) sum(x=='A',na.rm=T)/sum(!is.na(x)))
    atable=cbind(acnt, afrac)
    colnames(atable)=c('hmm.count.alive', 'hmm.frac.alive')
    rm(acnt, afrac)
    oligo.stats= data.frame(oligo.stats, atable[match(oligo.stats$GENEID, rownames(atable)),], stringsAsFactors=F)

    #-------------------------------------------------------------------
    plot.dir='/media/jbloom/d1/coupled_CRISPR/plots/HMM_v9/'
    makeHMMplots(plot.dir, oligo.stats, big.mm, conservation, dgsplit)
    
# ====================================================================================================

# confident alive oligos for example 
#    (oligo.stats[which(-log10(oligo.stats$ALL.p.value)>5& 
#                           oligo.stats$hmm2==1 & 
#                           oligo.stats$binarized.oligo.blup>0.5 &
#                           oligo.stats$hmm3==1)
#                                ,])
#    plot( oligo.stats$binarized.oligo.blup,  -log10(oligo.stats$ALL.p.value))
# Recalculate BLUP oligo and gene model and remove the essential genes 
# scale(guide.GCcontent)+scale(PAMvariantCNT)+
# BINARIZED
    formula.input=as.formula(slope.binarized~cnt+p_intercept+expt+binot+U6.Terminator+Score+guide.GCcontent+PAMvariantCNT+(1|GENEID)+(1|oligo))
    oligo.stats=doGLMER.oligo.gene(formula.input, data.in=big.mm[!big.mm$drop,], oligo.stats, lab='EssentialAll')

    formula.input=as.formula(slope.binarized~cnt+p_intercept+binot+U6.Terminator+Score+guide.GCcontent+PAMvariantCNT+(1|GENEID)+(1|oligo))
    oligo.stats=doGLMER.oligo.gene(formula.input, data.in=big.mm[!big.mm$drop & (big.mm$expt == 'WT'),], oligo.stats, lab='EssentialWT')

    formula.input=as.formula(slope.binarized~cnt+p_intercept+binot+U6.Terminator+Score+guide.GCcontent+PAMvariantCNT+(1|GENEID)+(1|oligo))
    oligo.stats=doGLMER.oligo.gene(formula.input, data.in=big.mm[!big.mm$drop & (big.mm$expt == 'NMD'),], oligo.stats, lab='EssentialNMD')

    formula.input=as.formula(slope.binarized~cnt+p_intercept+binot+U6.Terminator+Score+guide.GCcontent+PAMvariantCNT+(1|GENEID)+(1|oligo))
    oligo.stats=doGLMER.oligo.gene(formula.input, data.in=big.mm[!big.mm$drop & big.mm$domain.downstream,], oligo.stats, lab='Essential.DomainDownstream')

    formula.input=as.formula(slope.binarized~cnt+p_intercept+binot+U6.Terminator+Score+guide.GCcontent+PAMvariantCNT+(1|GENEID)+(1|oligo))
    oligo.stats=doGLMER.oligo.gene(formula.input, data.in=big.mm[!big.mm$drop & !big.mm$domain.downstream,], oligo.stats, lab='Essential.noDomainDownstream')

#slopes 
    formula.input=as.formula(slope~cnt+p_intercept+expt+binot+U6.Terminator+Score+guide.GCcontent+PAMvariantCNT+(1|GENEID)+(1|oligo))
    oligo.stats=doLMER.oligo.gene(formula.input, data.in=big.mm[!big.mm$drop,], oligo.stats, lab='EssentialAll')

    formula.input=as.formula(slope~cnt+p_intercept+binot+U6.Terminator+Score+guide.GCcontent+PAMvariantCNT+(1|GENEID)+(1|oligo))
    oligo.stats=doLMER.oligo.gene(formula.input, data.in=big.mm[!big.mm$drop & (big.mm$expt == 'WT'),], oligo.stats, lab='EssentialWT')

    formula.input=as.formula(slope~cnt+p_intercept+binot+U6.Terminator+Score+guide.GCcontent+PAMvariantCNT+(1|GENEID)+(1|oligo))
    oligo.stats=doLMER.oligo.gene(formula.input, data.in=big.mm[!big.mm$drop & (big.mm$expt == 'NMD'),], oligo.stats, lab='EssentialNMD')

    formula.input=as.formula(slope~cnt+p_intercept+binot+U6.Terminator+Score+guide.GCcontent+PAMvariantCNT+(1|GENEID)+(1|oligo))
    oligo.stats=doLMER.oligo.gene(formula.input, data.in=big.mm[!big.mm$drop & big.mm$domain.downstream,], oligo.stats, lab='Essential.DomainDownstream')

    formula.input=as.formula(slope~cnt+p_intercept+binot+U6.Terminator+Score+guide.GCcontent+PAMvariantCNT+(1|GENEID)+(1|oligo))
    oligo.stats=doLMER.oligo.gene(formula.input, data.in=big.mm[!big.mm$drop & !big.mm$domain.downstream,], oligo.stats, lab='Essential.noDomainDownstream')

    oligo.stats=oligo.stats[order(oligo.stats$unique.Index),]

#scale(guide.GCcontent)+scale(PAMvariantCNT)
#formula.input=as.formula(slope.binarized~cnt+p_intercept+expt+binot+U6.Terminator+Score+guide.GCcontent+PAMvariantCNT+(1|oligo))
#oligo.stats=doGLMER.oligo.gene(formula.input, data.in=big.mm[!big.mm$drop,], oligo.stats, lab='EssentialAll_noGene', doGene=FALSE)
#============================================================================================
save(oligo.stats, file=paste0(out.base.dir, 'processed/RData/oligo.stats.EssentialBlups.RData'))
load(paste0(out.base.dir, 'processed/RData/oligo.stats.EssentialBlups.RData'))



# BIG model with all----------------------------------------------------------------------------------- 
b3=big.mm[!big.mm$drop , ] # & (big.mm$expt == 'WT') , ]#& big.mm$dist_from_CDS_end<300,]
b3$expt=as.factor(b3$expt)
b3$oligo=as.factor(b3$oligo)
b3$domain.downstream=as.factor(b3$domain.downstream)
# scale(H0_w)+
#splic.genes=unlist(strsplit(GO.out[[1]][1,]$genes.in.set.sorted, ' '))
#spgf=b3$GENEID %in% splic.genes
s2=glmer(slope.binarized~cnt+expt+p_intercept+guide.GCcontent+dist_from_CDS_end+CDS_length+PAMvariantCNT+
                binot+U6.Terminator+Score+lcd+evolvability+mean.aa.perfect.conserved+viable_annotation+
                end.conservation+domain.downstream+(1|GENEID)+(1|oligo), data=b3, verbose=T,  family=binomial(logit), 
                ,control = glmerControl(optCtrl= list(maxfun=5e5), calc.derivs=FALSE ))
#save(s2, file=paste0(out.base.dir, 'processed/RData/bigModel.RData'))
load(paste0(out.base.dir, 'processed/RData/bigModel.RData'))
s2a=Anova(s2, type='III') 
summary(s2)
s2a=s2a[order(s2a$Chisq, decreasing=T),]
WriteXLS(s2a, '/home/jbloom/Dropbox/Public/CoupledCRISPR/bigRegressionBinarized.xls', row.names=T )


### GO enrichment
testGenes=oligo.stats$binarized.gene.blup.EssentialAll[match(unique(oligo.stats[!oligo.stats$drop,]$GENEID), oligo.stats$GENEID)]
names(testGenes)=oligo.stats$GENEID[match(unique(oligo.stats[!oligo.stats$drop,]$GENEID), oligo.stats$GENEID)]
testGenes=sort(testGenes, decreasing=T)
GO.out=doGO(testGenes)
WriteXLS(GO.out, '/home/jbloom/Dropbox/Public/CoupledCRISPR/GO_gene_blups_ks_test.xls', SheetNames=names(GO.out))

o10=oligo.stats[!oligo.stats$drop,]
#gset=head(unique(o10$GENEID[order(o10$binarized.gene.blup, decreasing=T)]), 15) 
#osub=oligo.stats[oligo.stats$hmm2=='A' &!oligo.stats$drop , ] # & oligo.stats$binarized.oligo.blup>0 & !oligo.stats$dubious, ] 
oss=split(o10, o10$GENEID)
fracA=sapply(oss, function(x) x$hmm.frac.alive[1])
fracA=sort(fracA, decreasing=T)
GO.outHMM=doGO(fracA)
WriteXLS(GO.outHMM, '/home/jbloom/Dropbox/Public/CoupledCRISPR/GO_fracAlive_ks_test.xls', SheetNames=names(GO.out))


#frac1=sapply(oss, function(x) sum(x$hmm2=='A',na.rm=T)>0)
#frac1GO.E=doGO.threshold(names(which(frac1)), names(frac1))
#frac1GO.D=doGO.threshold(names(which(!frac1)), names(frac1))
#WriteXLS(frac1GO.E, '/home/jbloom/Desktop/LabMeeting-041217/GO_1e.xls', SheetNames=names(frac1GO.E))
#WriteXLS(frac1GO.D, '/home/jbloom/Desktop/LabMeeting-041217/GO_1d.xls', SheetNames=names(frac1GO.D))
#frac2=sapply(oss, function(x) sum(x$hmm2=='A',na.rm=T)>2)
#frac5=sapply(oss, function(x) sum(x$hmm2=='A',na.rm=T)>4)
#frac5GO.E=doGO.threshold(names(which(frac5)), names(frac5))
#---------------------------------------------------------------------------------------------------------


on2=oligo.stats[match(unique(oligo.stats$GENEID), oligo.stats$GENEID),]


# Figure 3 -----------------------------------------------------------------------------------------------
column="binarized.oligo.blup.EssentialAll"
column="binarized.oligo.blup.Essential.DomainDownstream"
column="binarized.oligo.blup.Essential.noDomainDownstream"
column="binarized.oligo.blup.EssentialWT"
column="binarized.oligo.blup.EssentialNMD"
pdf(file='/home/jbloom/Dropbox/Public/CoupledCRISPR/Figures/Figure3_150.pdf', width=10, height=10)
pdf(file='/home/jbloom/Dropbox/Public/CoupledCRISPR/Figures/Figure3_150_WTvNMD.pdf', width=10, height=10)
pdf(file='/home/jbloom/Dropbox/Public/CoupledCRISPR/Figures/Figure3_150_DomaninvNoDomain.pdf', width=10, height=10)
par(mfrow=c(2,1))
plot(oligo.stats$dist_from_CDS_end, oligo.stats[,column], xlim=c(150,0),xaxt='n', ylab='PTC tolerance', xlab='AA distance from C-terminal', col='#00000044', main=column, cex=.75)
axis(1, at=rev(seq(0,150,5)))
x2=na.omit(data.frame(distx=oligo.stats$dist_from_CDS_end, y=oligo.stats[,column]))
x3=x2[x2$distx<500,]
s3=segmented(lm(y~distx, data=x3), seg.Z=~distx, control=seg.control(n.boot=100))
seg.fit=cbind(x3$distx, predict(s3))
seg.fit=seg.fit[order(seg.fit[,1]),]
points(seg.fit[,1], seg.fit[,2], col='blue', type='l', lwd=4)
abline(v=confint.segmented(s3)$distx, col='blue')
#
ods=split(oligo.stats[,column], oligo.stats$dist_from_CDS_end)
odsq=t(sapply(ods, function(x) quantile(x, c(.25,.5,.75), na.rm=T)))
no=as.numeric(rownames(odsq))
segments(no-.5, odsq[,1], no+.5, odsq[,1], lwd=2)
segments(no-.5, odsq[,2], no+.5, odsq[,2], col='red', lwd=2)
segments(no-.5, odsq[,3], no+.5, odsq[,3], lwd=2)
     
dev.off()
#---------------------------------------------------------------------------------------------------------





# Figure 1 pilot testing ------------------------------------------------------------------------------------
    git.dir='/media/jbloom/d1/coupled_CRISPR/Experiments/Preliminary/Pilot_092015/'
    load(paste0(git.dir, 'summaryStats.RData'))
    #load(paste0(git.dir, 'sampleTable.RData'))

    #summary.stats.table=t(do.call('rbind', summary.stats))
    summary.stats.table=summary.stats
    pilot.norm=(summary.stats.table)/rowSums(summary.stats.table)
    pilot.norm[,1]=pilot.norm[,1]+pilot.norm[,'X32']
    pilot.norm=pilot.norm[,-match('X32', colnames(pilot.norm))]
    pilot.norm[,'X31']=pilot.norm[,'X31']+pilot.norm[,'X33']
    pilot.norm=pilot.norm[,-match('X33', colnames(pilot.norm))]

    #1 = matches haplotype no mismatches or indels in 50bp window of pam
    #2 = indel in 50 bp window around pam
    #4 = WT
    #31 = partial edit (+ mismatch)
    #32 = full edit (+ mismatch)
    #33 = no edit (+ mismatch)

    pilot.norm.drop=t(pilot.norm[rev(c(1:6,8, 27)),])
    pdf(file='/home/jbloom/Dropbox/Public/CoupledCRISPR/Figures/Figure1B.pdf', width=22, height=8)
    par(oma=c(1,20,1,1))
    bp=barplot(pilot.norm.drop, col=c('black','red','white','grey'), # density=c(100,55,25,0),
               horiz=TRUE, xlim=c(0,1), yaxt='n', 
               beside=T, xlab='fraction of reads', space=c(0,3), width=2 )#, main='Figure1B')
    axis(2, bp[2,], colnames(pilot.norm.drop), las=2)
    legend('bottomright',fill=c('black','red',  'white','grey'), #density= #rev(c(100,55,25,0)) , 
                                               (c('matches expected edit', 
                                                    'indel near edit',
                                                    'unedited', 
                                                    'partial edit')))
    dev.off()

    # supplement 1B (contrast NEJ1delete vs WT)
    pilot.norm.drop=t(pilot.norm[rev(c(20:23)),])
    pdf(file='/home/jbloom/Dropbox/Public/CoupledCRISPR/Figures/SFigure1B.pdf', width=22, height=8)
    par(oma=c(1,20,1,1))
    bp=barplot(pilot.norm.drop, col=c('black','red', 'white','grey'), # density=c(100,55,25,0),
               horiz=TRUE, xlim=c(0,1), yaxt='n', 
               beside=T, xlab='fraction of reads', space=c(0,3), width=2 )#, main='Figure1B')
    axis(2, bp[2,], colnames(pilot.norm.drop), las=2)
    legend('bottomright',fill=c('black','red',  'white', 'grey'), #density= #rev(c(100,55,25,0)) , 
                                               (c('matches expected edit', 
                                                    'indel near edit',
                                                    'unedited', 
                                                    'partial edit')))
    dev.off()
    # NMD delete
    pilot.norm.drop=t(pilot.norm[rev(c(9:14,16)),])
    pdf(file='/home/jbloom/Dropbox/Public/CoupledCRISPR/Figures/SFigure1B_NMDm.pdf', width=22, height=8)
    par(oma=c(1,20,1,1))
    bp=barplot(pilot.norm.drop, col=c('black','red', 'white','grey'), # density=c(100,55,25,0),
               horiz=TRUE, xlim=c(0,1), yaxt='n', 
               beside=T, xlab='fraction of reads', space=c(0,3), width=2 )#, main='Figure1B')
    axis(2, bp[2,], colnames(pilot.norm.drop), las=2)
    legend('bottomright',fill=c('black','red',  'white','grey'), #density= #rev(c(100,55,25,0)) , 
                                               (c('matches expected edit', 
                                                    'indel near edit',
                                                    'unedited', 
                                                    'partial edit')))
    dev.off()
#-------------------------------------------------------------------------------------------------------------

# Figure 2 slopes (Dubious vs essential)---------------------------------------------------------------------
    pdf(file='/home/jbloom/Dropbox/Public/CoupledCRISPR/Figures/Figure2.pdf', width=12, height=12)
    par(oma=c(1,1,1,1))
    test= stripchart(oligo.stats$binarized.oligo.blup.ng~oligo.stats$dubious, vertical=T, method='jitter', col='#00000012', pch=19,
               group.names=c('essential', 'dubious'),
               ylab='PTC tolerance', xlim=c(0.5,2.5) )
    boxplot(oligo.stats$binarized.oligo.blup.ng~oligo.stats$dubious, add=T, xlab='', names='', ylab='', outline=F) 
    dev.off()

    dup.guides=oligo.stats$guide[which(duplicated(oligo.stats$guide))]
    idg= oligo.stats$guide %in% dup.guides
    osgs=split(oligo.stats[idg,], oligo.stats[idg,]$guide)
    osgs.table=t(sapply(osgs, function(x) c(x$binarized.oligo.blup.ng[x$dubious], x$binarized.oligo.blup.ng[!x$dubious])))
    pdf(file='/home/jbloom/Dropbox/Public/CoupledCRISPR/Figures/sFigure2_matched_gRNA.pdf', width=8, height=8)
    plot(osgs.table, xlim=c(-2,4), ylim=c(-2,4), ylab='essential PTC tolerance', xlab='dubious PTC tolerance' )
    abline(0,1)
    dev.off()
    t.test(osgs.table[,1], osgs.table[,2], paired=T)

#------------------------------------------------------------------------------------------------------------































#gAp=plot(gammALL$gam)
#gApD=plot(gammALL.dd$gam)
#gApnD=plot(gammALL.ndd$gam)
#plot(gAp[[1]]$x, gAp[[1]]$fit, type='l', xlim=c(200,0))
#points(gAp[[1]]$x, gAp[[1]]$fit+gAp[[1]]$se, type='l', lty=2)
#points(gAp[[1]]$x, gAp[[1]]$fit-gAp[[1]]$se, type='l', lty=2)
#om=oligo.stats[match(unique(oligo.stats$GENEID), oligo.stats$GENEID),]
#hist(om$binarized.gene.blup)
#length(unique(oligo.stats$GENEID[oligo.stats$binarized.oligo.blup>0]))

## Figure 3 distance from end  (WT vs NMD) potentially also (DOMAIM vs no DOMAIN)---------------------------------------------------------------
    xlims=150
    pdf(file='~/Desktop/Figure3_150.pdf', width=10, height=10)
    x=oligo.stats[!oligo.stats$drop & !is.na(oligo.stats$binarized.oligo.blup.EssentialAll) & oligo.stats$dist_from_CDS_end<xlims,]    
    par(mfrow=c(2,1),oma=c(1,2,1,1))
    plot(x$dist_from_CDS_end, x$binarized.oligo.blup.EssentialWT, xlim=c(xlims,0), col='#0000ff33', xlab="amino acids from 3' end of protein", ylab="log[(odds PTC surviving)/(odds PTC dying)]" )#, main='Figure3')
    points(x$dist_from_CDS_end, x$binarized.oligo.blup.EssentialNMD, xlim=c(xlims,0), col='#00800033', xlab="amino acids from 3' end of protein" )
    xx=cbind(x$dist_from_CDS_end, x$binarized.oligo.blup.EssentialWT)
    xx=na.omit(xx)
    data=xx
    sc=spline.cis(data,B=1000, alpha=.01)
    lines(x=sc$x, y=sc$main.curve, col='blue', lwd=2)
    lines(x=sc$x, y=sc$lower.ci, col='blue', lty=2)
    lines(x=sc$x, y=sc$upper.ci, col='blue', lty=2)

    xx=cbind(x$dist_from_CDS_end, x$binarized.oligo.blup.EssentialNMD)
    xx=na.omit(xx)
    data=xx
    sc=spline.cis(data,B=1000, alpha=.01)
    lines(x=sc$x, y=sc$main.curve, col='green', lwd=2)
    lines(x=sc$x, y=sc$lower.ci, col='green', lty=2)
    lines(x=sc$x, y=sc$upper.ci, col='green', lty=2)
    legend('topleft', fill=c('blue', 'green'), legend=c('WT', 'NMDm'))


    plot(x$dist_from_CDS_end[x$domain.downstream], x$binarized.oligo.blup.EssentialAll[x$domain.downstream], xlim=c(xlims,0), col='#00000033', xlab="amino acids from 3' end of protein", 
         ylab="log[(odds PTC surviving)/(odds PTC dying)]" ) #, main='Figure3')
    points(x$dist_from_CDS_end[!x$domain.downstream], x$binarized.oligo.blup.EssentialAll[!x$domain.downstream], xlim=c(xlims,0), col='#ff000033', xlab="amino acids from 3' end of protein" )
    xx=cbind(x$dist_from_CDS_end[x$domain.downstream], x$binarized.oligo.blup.EssentialAll[x$domain.downstream])
    xx=na.omit(xx)
    data=xx
    sc=spline.cis(data,B=1000, alpha=.01)
    lines(x=sc$x, y=sc$main.curve, col='black', lwd=2)
    lines(x=sc$x, y=sc$lower.ci, col='black', lty=2)
    lines(x=sc$x, y=sc$upper.ci, col='black', lty=2)
    xx=cbind(x$dist_from_CDS_end[!x$domain.downstream], x$binarized.oligo.blup.EssentialAll[!x$domain.downstream])
    xx=na.omit(xx)
    data=xx
    sc=spline.cis(data,B=1000, alpha=.01)
    lines(x=sc$x, y=sc$main.curve, col='red', lwd=2)
    lines(x=sc$x, y=sc$lower.ci, col='red', lty=2)
    lines(x=sc$x, y=sc$upper.ci, col='red', lty=2)
    legend('topleft', fill=c('black', 'red'), legend=c('PTCs with domain downstream', 'PTCs with no domain downstream'))
    title(outer=T , 'Figure3')
    dev.off()

#----------------------------------------------------------------------------------------------------------------------------------------------
# HMM based analysis 
e2=cbind(oligo.stats$hmm2, oligo.stats$dist_from_CDS_end)#[!oligo.stats$dubious,]
e2[,1]=as.numeric(e2[,1])-1
e2=na.omit(e2)
se2=split(e2[,1], e2[,2])

se2.frac=sapply(se2, function(x) sum(x)/length(x))
plot(as.numeric(names(se2.frac)), se2.frac, xlim=c(1000,0), cex=sqrt(sapply(se2, length)/10))

smb=matrix(0, length(se2.frac), 100)
for(i in 1:100) {
smb[,i]=sapply(se2, function(x) {y=sample(x, replace=T);  return(sum(y)/length(y)) })
}
segments(as.numeric(names(se2.frac)),
         apply(smb,1, quantile, .05),
         as.numeric(names(se2.frac)),
         apply(smb,1, quantile, .95), col='red')

##### data subset without domain downstream
par(mfrow=c(2,1))
e3=cbind(oligo.stats$hmm2, oligo.stats$dist_from_CDS_end)[!oligo.stats$domain.downstream,]
e3[,1]=as.numeric(e3[,1])-1
e3=na.omit(e3)
se3=split(e3[,1], e3[,2])
se3.frac=sapply(se3, function(x) sum(x)/length(x))
plot(as.numeric(names(se3.frac)), se3.frac, xlim=c(150,0), cex=sqrt(sapply(se2, length)/20), 
     main='no domain downstream',
     xlab="PTC aa position relative to 3'",
     ylab='fraction of PTCs tolerated')

smb3=matrix(0, length(se3.frac), 100)
for(i in 1:100) {
smb3[,i]=sapply(se3, function(x) {y=sample(x, replace=T);  return(sum(y)/length(y)) })
}
segments(as.numeric(names(se3.frac)),
         apply(smb3,1, quantile, .05),
         as.numeric(names(se3.frac)),
         apply(smb3,1, quantile, .95), col='red')

##### data subset with domain downstream
e4=cbind(oligo.stats$hmm2, oligo.stats$dist_from_CDS_end)[oligo.stats$domain.downstream,]
e4[,1]=as.numeric(e4[,1])-1
e4=na.omit(e4)
se4=split(e4[,1], e4[,2])
se4.frac=sapply(se4, function(x) sum(x)/length(x))

plot(as.numeric(names(se4.frac)), se4.frac, xlim=c(150,0), cex=sqrt(sapply(se2, length)/20), 
     main='domain downstream',
     xlab="PTC aa position relative to 3'",
     ylab='fraction of PTCs tolerated')

smb4=matrix(0, length(se4.frac), 100)
for(i in 1:100) {
smb4[,i]=sapply(se4, function(x) {y=sample(x, replace=T);  return(sum(y)/length(y)) })
}
segments(as.numeric(names(se4.frac)),
         apply(smb4,1, quantile, .05),
         as.numeric(names(se4.frac)),
         apply(smb4,1, quantile, .95), col='red')




frac6=sapply(oss, function(x) sum(x$hmm2=='A',na.rm=T)>5)
frac7=sapply(oss, function(x) sum(x$hmm2=='A',na.rm=T)>6)
points(as.numeric(names(se2.frac)), apply(smb,1, quantile, .95), col='red')
points(as.numeric(names(se2.frac)), apply(smb,1, quantile, .05), col='red')
apply(smb,1, quantile, .5)
apply(smb,1, quantile, .1)


#bs.se2=matrreplicate(100, function(i) { sapply(se2, function(x) {y=sample(x);  return(sum(y)/length(y)) })  })
#, xlim=c(150,0))

plot(e2[,2], predict(glm(e2[,1]~e2[,2], family=binomial(link='logit'))))
plot(e2[,2], predict(glm(e2[,1]~e2[,2], family=binomial(link='logit')), type='response'))
        
        
# fit binarized slope model  Figure or table 4
# ALL DATA--------------------------------------------------------------
    png(file= '/home/jbloom/Desktop/LabMeeting-041217/02_ranef_slope_models.png', width=1080, height=1080)
    b3=big.mm[!big.mm$drop , ] # & (big.mm$expt == 'WT') , ]#& big.mm$dist_from_CDS_end<300,]
    b3$expt=as.factor(b3$expt)
    b3$domain.downstream=as.factor(b3$domain.downstream)
    b3$oligo=as.factor(b3$oligo)
    # scale(H0_w)+
    #scale(CDS_length)+
    #binot+
    # subset of data out to 125
    # 2 step ... fit full model... then remove terms that explain nothin 
    #mrna or some other data set
    #splic.genes=unlist(strsplit(GO.out[[1]][1,]$genes.in.set.sorted, ' '))
    #spgf=b3$GENEID %in% splic.genes
    
    s2=lmer(slope~ cnt+
                    as.factor(expt)+
                    p_intercept+
                    guide.GCcontent+
                    #dist_from_CDS_end+
                    CDS_length+
                    PAMvariantCNT+
                    binot+U6.Terminator+Score+
                    #lcd+
                    #evolvability+
                    #mean.aa.perfect.conserved+
                    #viable_annotation+
                    #end.conservation+
                    #domain.downstream+
                    #spgf+
                    (1|GENEID)+(1|oligo), data=b3, verbose=T)
    s2a=(Anova(s2, type='III'))
    s2resid=residuals(s2)
    s2ranef=ranef(s2)
    #plot(jitter(b3$dist_from_CDS_end),b3$slope, xlim=c(150,0), col="#00000044",ylim=c(-.125, .05),
    #     xlab=c("AA distance from 3' end"), ylab='barcode slope')
    plot(oligo.stats$dist_from_CDS_end[match(rownames(s2ranef$oligo), oligo.stats$oligo)], s2ranef$oligo[,1], xlim=c(100,0), col="#000000BB",
            xlab=c("AA distance from 3' end"), ylab='oligo random effect', ylim=c(-.04,.04))
    x2=cbind(oligo.stats$dist_from_CDS_end[match(rownames(s2ranef$oligo), oligo.stats$oligo)], s2ranef$oligo[,1])
    colnames(x2)=c('distx', 'y')
    x2=data.frame(x2)
    #for(s in seq(10,1000,10)){
    x3=x2[x2$distx<500,]
    l1=lm(y~distx, data=x3)
    s3=segmented(l1, seg.Z=~distx, control=seg.control(n.boot=100))
    s3n=nls(y~c+a*(exp(-b*distx)), data=x3, start=list(c=-.001,a=0.001, b=.05))
    which.min(c(BIC(s3n), BIC(s3)))
    seg.fit=cbind(x3$distx, predict(s3))
    seg.fit=seg.fit[order(seg.fit[,1]),]

    seg.fitn=cbind(x3$distx, predict(s3n))
    seg.fitn=seg.fitn[order(seg.fitn[,1]),]

    points(seg.fit[,1], seg.fit[,2], col='blue', type='l', lwd=3)
    points(seg.fit[,1], seg.fitn[,2], col='red', type='l', lwd=3)
    abline(v=confint.segmented(s3)$distx, col='blue')

    text(90, .04, expression(y = -0.00129+0.01469(e^(-0.05920*x))), col='red', cex=1.5)
    dev.off()
    xx=1:150
    yy=-0.00129+0.01469*(exp(-0.05920*xx))

#---------------------------------------------------------------------
# domain downstream
    png(file= '/home/jbloom/Desktop/LabMeeting-041217/02_ranef_slope_models+by_domain.png', width=1080, height=1080)
    par(mfrow=c(2,1))
    b3=big.mm[!big.mm$drop & big.mm$domain.downstream, ] # & (big.mm$expt == 'WT') , ]#& big.mm$dist_from_CDS_end<300,]
    b3$expt=as.factor(b3$expt)
    b3$domain.downstream=as.factor(b3$domain.downstream)
    b3$oligo=as.factor(b3$oligo)
    # scale(H0_w)+
    #scale(CDS_length)+
    #binot+
    # subset of data out to 125
    # 2 step ... fit full model... then remove terms that explain nothin 
    #mrna or some other data set
    #splic.genes=unlist(strsplit(GO.out[[1]][1,]$genes.in.set.sorted, ' '))
    #spgf=b3$GENEID %in% splic.genes
    s2=lmer(slope~ cnt+
                    as.factor(expt)+
                    p_intercept+
                    guide.GCcontent+
                    #dist_from_CDS_end+
                    CDS_length+
                    PAMvariantCNT+
                    binot+U6.Terminator+Score+
                    #lcd+
                    #evolvability+
                    #mean.aa.perfect.conserved+
                    #viable_annotation+
                    #end.conservation+
                    #domain.downstream+
                    #spgf+
                    (1|GENEID)+(1|oligo), data=b3, verbose=T)
    s2a=(Anova(s2, type='III'))
    s2resid=residuals(s2)
    s2ranef=ranef(s2)
    #plot(jitter(b3$dist_from_CDS_end),b3$slope, xlim=c(150,0), col="#00000044",ylim=c(-.125, .05),
    #     xlab=c("AA distance from 3' end"), ylab='barcode slope')
    plot(oligo.stats$dist_from_CDS_end[match(rownames(s2ranef$oligo), oligo.stats$oligo)], s2ranef$oligo[,1], xlim=c(100,0), col="#00000088",
            xlab=c("AA distance from 3' end"), ylab='oligo random effect', ylim=c(-.04,.04), main='has downtream domain')
    x2=cbind(oligo.stats$dist_from_CDS_end[match(rownames(s2ranef$oligo), oligo.stats$oligo)], s2ranef$oligo[,1])
    colnames(x2)=c('distx', 'y')
    x2=data.frame(x2)
    #for(s in seq(10,1000,10)){
    x3=x2[x2$distx<500,]
    l1=lm(y~distx, data=x3)
    s3=segmented(l1, seg.Z=~distx, control=seg.control(n.boot=100))
    s3n=nls(y~c+a*(exp(-b*distx)), data=x3, start=list(c=-.001,a=0.001, b=.05))
    which.min(c(BIC(s3n), BIC(s3)))
    seg.fit=cbind(x3$distx, predict(s3))
    seg.fit=seg.fit[order(seg.fit[,1]),]

    seg.fitn=cbind(x3$distx, predict(s3n))
    seg.fitn=seg.fitn[order(seg.fitn[,1]),]

    points(seg.fit[,1], seg.fit[,2], col='blue', type='l', lwd=3)
    points(seg.fit[,1], seg.fitn[,2], col='red', type='l', lwd=3)
    #abline(v=confint.segmented(s3)$distx, col='blue')

    #### no domain downstream
b3=big.mm[!big.mm$drop & !big.mm$domain.downstream, ] # & (big.mm$expt == 'WT') , ]#& big.mm$dist_from_CDS_end<300,]
    b3$expt=as.factor(b3$expt)
    b3$domain.downstream=as.factor(b3$domain.downstream)
    b3$oligo=as.factor(b3$oligo)
    # scale(H0_w)+
    #scale(CDS_length)+
    #binot+
    # subset of data out to 125
    # 2 step ... fit full model... then remove terms that explain nothin 
    #mrna or some other data set
       s2=lmer(slope~ cnt+
                    as.factor(expt)+
                    p_intercept+
                    guide.GCcontent+
                    #dist_from_CDS_end+
                    CDS_length+
                    PAMvariantCNT+
                    binot+U6.Terminator+Score+
                    #lcd+
                    #evolvability+
                    #mean.aa.perfect.conserved+
                    #viable_annotation+
                    #end.conservation+
                    #domain.downstream+
                    #spgf+
                    (1|GENEID)+(1|oligo), data=b3, verbose=T)
    s2a=(Anova(s2, type='III'))
    s2resid=residuals(s2)
    s2ranef=ranef(s2)
    #plot(jitter(b3$dist_from_CDS_end),b3$slope, xlim=c(150,0), col="#00000044",ylim=c(-.125, .05),
    #     xlab=c("AA distance from 3' end"), ylab='barcode slope')
    #x11()
    plot(oligo.stats$dist_from_CDS_end[match(rownames(s2ranef$oligo), oligo.stats$oligo)], s2ranef$oligo[,1], xlim=c(100,0), col="#00000088",
            xlab=c("AA distance from 3' end"), ylab='oligo random effect', ylim=c(-.04,.04), main='no downtream domain')
    x2=cbind(oligo.stats$dist_from_CDS_end[match(rownames(s2ranef$oligo), oligo.stats$oligo)], s2ranef$oligo[,1])
    colnames(x2)=c('distx', 'y')
    x2=data.frame(x2)
    #for(s in seq(10,1000,10)){
    x3=x2[x2$distx<500,]
    l1=lm(y~distx, data=x3)
    s3=segmented(l1, seg.Z=~distx, control=seg.control(n.boot=100))
    s3n=nls(y~c+a*(exp(-b*distx)), data=x3, start=list(c=-.001,a=0.001, b=.05))
    which.min(c(BIC(s3n), BIC(s3)))
    seg.fit=cbind(x3$distx, predict(s3))
    seg.fit=seg.fit[order(seg.fit[,1]),]

    seg.fitn=cbind(x3$distx, predict(s3n))
    seg.fitn=seg.fitn[order(seg.fitn[,1]),]
    points(seg.fit[,1], seg.fit[,2], col='blue', type='l', lwd=3)
    points(seg.fit[,1], seg.fitn[,2], col='red', type='l', lwd=3)
    abline(v=confint.segmented(s3)$distx, col='blue')
dev.off()

####BY WT and NMD
    png(file= '/home/jbloom/Desktop/LabMeeting-041217/02_ranef_slope_models+by_WT_vs_NMD.png', width=1080, height=1080)
    par(mfrow=c(2,1))
    b3=big.mm[!big.mm$drop & big.mm$expt == 'WT', ]#& big.mm$dist_from_CDS_end<300,]
    b3$expt=as.factor(b3$expt)
    b3$domain.downstream=as.factor(b3$domain.downstream)
    b3$oligo=as.factor(b3$oligo)
    # scale(H0_w)+
    #scale(CDS_length)+
    #binot+
    # subset of data out to 125
    # 2 step ... fit full model... then remove terms that explain nothin 
    #mrna or some other data set
    #splic.genes=unlist(strsplit(GO.out[[1]][1,]$genes.in.set.sorted, ' '))
    #spgf=b3$GENEID %in% splic.genes
    s2=lmer(slope~ cnt+
                    p_intercept+
                    guide.GCcontent+
                    #dist_from_CDS_end+
                    CDS_length+
                    PAMvariantCNT+
                    binot+U6.Terminator+Score+
                    #lcd+
                    #evolvability+
                    #mean.aa.perfect.conserved+
                    #viable_annotation+
                    #end.conservation+
                    #domain.downstream+
                    #spgf+
                    (1|GENEID)+(1|oligo), data=b3, verbose=T)
    s2a=(Anova(s2, type='III'))
    s2resid=residuals(s2)
    s2ranef=ranef(s2)
    #plot(jitter(b3$dist_from_CDS_end),b3$slope, xlim=c(150,0), col="#00000044",ylim=c(-.125, .05),
    #     xlab=c("AA distance from 3' end"), ylab='barcode slope')
    plot(oligo.stats$dist_from_CDS_end[match(rownames(s2ranef$oligo), oligo.stats$oligo)], s2ranef$oligo[,1], xlim=c(100,0), col="#00000088",
            xlab=c("AA distance from 3' end"), ylab='oligo random effect', ylim=c(-.04,.04), main='WT')
    x2=cbind(oligo.stats$dist_from_CDS_end[match(rownames(s2ranef$oligo), oligo.stats$oligo)], s2ranef$oligo[,1])
    colnames(x2)=c('distx', 'y')
    x2=data.frame(x2)
    #for(s in seq(10,1000,10)){
    x3=x2[x2$distx<500,]
    l1=lm(y~distx, data=x3)
    s3=segmented(l1, seg.Z=~distx, control=seg.control(n.boot=100))
    print(summary(s3))
    s3n=nls(y~c+a*(exp(-b*distx)), data=x3, start=list(c=-.001,a=0.001, b=.05))
    which.min(c(BIC(s3n), BIC(s3)))
    seg.fit=cbind(x3$distx, predict(s3))
    seg.fit=seg.fit[order(seg.fit[,1]),]

    seg.fitn=cbind(x3$distx, predict(s3n))
    seg.fitn=seg.fitn[order(seg.fitn[,1]),]

    points(seg.fit[,1], seg.fit[,2], col='blue', type='l', lwd=3)
    points(seg.fit[,1], seg.fitn[,2], col='red', type='l', lwd=3)
    abline(v=confint.segmented(s3)$distx, col='blue')

    b3=big.mm[!big.mm$drop & big.mm$expt == 'NMD', ]#& big.mm$dist_from_CDS_end<300,]
    b3$expt=as.factor(b3$expt)
    b3$domain.downstream=as.factor(b3$domain.downstream)
    b3$oligo=as.factor(b3$oligo)
    # scale(H0_w)+
    #scale(CDS_length)+
    #binot+
    # subset of data out to 125
    # 2 step ... fit full model... then remove terms that explain nothin 
    #mrna or some other data set
    #splic.genes=unlist(strsplit(GO.out[[1]][1,]$genes.in.set.sorted, ' '))
    #spgf=b3$GENEID %in% splic.genes
    s2=lmer(slope~ cnt+
                    p_intercept+
                    guide.GCcontent+
                    #dist_from_CDS_end+
                    CDS_length+
                    PAMvariantCNT+
                    binot+U6.Terminator+Score+
                    #lcd+
                    #evolvability+
                    #mean.aa.perfect.conserved+
                    #viable_annotation+
                    #end.conservation+
                    #domain.downstream+
                    #spgf+
                    (1|GENEID)+(1|oligo), data=b3, verbose=T)
    s2a=(Anova(s2, type='III'))
    s2resid=residuals(s2)
    s2ranef=ranef(s2)
    #plot(jitter(b3$dist_from_CDS_end),b3$slope, xlim=c(150,0), col="#00000044",ylim=c(-.125, .05),
    #     xlab=c("AA distance from 3' end"), ylab='barcode slope')
    plot(oligo.stats$dist_from_CDS_end[match(rownames(s2ranef$oligo), oligo.stats$oligo)], s2ranef$oligo[,1], xlim=c(100,0), col="#00000088",
            xlab=c("AA distance from 3' end"), ylab='oligo random effect', ylim=c(-.04,.04), main='NMD del')
    x2=cbind(oligo.stats$dist_from_CDS_end[match(rownames(s2ranef$oligo), oligo.stats$oligo)], s2ranef$oligo[,1])
    colnames(x2)=c('distx', 'y')
    x2=data.frame(x2)
    #for(s in seq(10,1000,10)){
    x3=x2[x2$distx<500,]
    l1=lm(y~distx, data=x3)
    s3=segmented(l1, seg.Z=~distx, control=seg.control(n.boot=100))
    print(summary(s3))

    s3n=nls(y~c+a*(exp(-b*distx)), data=x3, start=list(c=-.001,a=0.001, b=.05))
    which.min(c(BIC(s3n), BIC(s3)))
    seg.fit=cbind(x3$distx, predict(s3))
    seg.fit=seg.fit[order(seg.fit[,1]),]

    seg.fitn=cbind(x3$distx, predict(s3n))
    seg.fitn=seg.fitn[order(seg.fitn[,1]),]

    points(seg.fit[,1], seg.fit[,2], col='blue', type='l', lwd=3)
    points(seg.fit[,1], seg.fitn[,2], col='red', type='l', lwd=3)
    abline(v=confint.segmented(s3)$distx, col='blue')
    dev.off()









### ADD Back terms 




b3=big.mm[!big.mm$drop , ] # & (big.mm$expt == 'WT') , ]#& big.mm$dist_from_CDS_end<300,]
b3$expt=as.factor(b3$expt)
#b3$domain.downstream=as.factor(b3$domain.downstream)
b3$oligo=as.factor(b3$oligo)
# scale(H0_w)+
#scale(CDS_length)+
#binot+
# subset of data out to 125
# 2 step ... fit full model... then remove terms that explain nothin 
#mrna or some other data set
splic.genes=unlist(strsplit(GO.out[[1]][1,]$genes.in.set.sorted, ' '))
spgf=b3$GENEID %in% splic.genes
s2=lmer(slope~cnt+
                as.factor(expt)+
                p_intercept+
                guide.GCcontent+
                dist_from_CDS_end+
                CDS_length+
                PAMvariantCNT+
                binot+U6.Terminator+Score+
                lcd+
                evolvability+
                mean.aa.perfect.conserved+
                viable_annotation+
                end.conservation+
                domain.downstream+
                spgf+
                #spgf+
                (1|GENEID)+(1|oligo), data=b3, verbose=T)
s2a2=anova(s2)
s2a=Anova(s2, type='III') #, test.statistic=c("F"))
s2a[order(s2a$Chisq, decreasing=T),]


# DO GO 




s2resid=residuals(s2)
s2ranef=ranef(s2)

plot(jitter(b3$dist_from_CDS_end),b3$slope, xlim=c(150,0), col="#00000044",ylim=c(-.125, .05),
     xlab=c("AA distance from 3' end"), ylab='barcode slope')
plot(oligo.stats$dist_from_CDS_end[match(rownames(s2ranef$oligo), oligo.stats$oligo)], s2ranef$oligo[,1], xlim=c(150,0), col="#00000055")
x2=cbind(oligo.stats$dist_from_CDS_end[match(rownames(s2ranef$oligo), oligo.stats$oligo)], s2ranef$oligo[,1])
colnames(x2)=c('distx', 'y')
#x2=data.frame(x2)
#for(s in seq(10,1000,10)){
x3=x2[x2$distx<500,]
l1=lm(y~distx, data=x3)
s3=segmented(l1, seg.Z=~distx)
lines(s3, col='red')






















plot_by_pos=with(oligo.stats[!oligo.stats$drop & !is.na(oligo.stats$binarized.oligo.blup.EssentialAll) & oligo.stats$dist_from_CDS_end<150,], {
   

            
    #blist=split(x$binarized.oligo.blup.EssentialAll[x$domain.downstream],x$dist_from_CDS_end[x$domain.downstream])
    #blist2=split(x$binarized.oligo.blup.EssentialAll[!x$domain.downstream],x$dist_from_CDS_end[!x$domain.downstream])



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
#--------------------------------------------------------------------------------------------------------------





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
#hist(oligo.stats$binarized.gene.blup[oligo.stats$dubious])
obthresh=quantile(oligo.stats$binarized.gene.blup[oligo.stats$dubious], .5)
geneschosen=unique(oligo.stats$GENEID[oligo.stats$binarized.gene.blup>obthresh & !oligo.stats$drop ])


# constrain background set
GO.out.threshold=doGO.threshold(geneschosen, unique(oligo.stats$GENEID[!oligo.stats$dubious]))
fdr.out=qvalue(unlist(sapply(GO.out.threshold, function(x) x$p.value)),fdr.level=.01)  #lambda=seq(.3,.5,.05), fdr.level=.01)
WriteXLS(GO.out.threshold, '/home/jbloom/Dropbox/Public/CoupledCRISPR/GO_gene_blups_binary_threshold.xls', SheetNames=names(GO.out.threshold))

# end of gene conservation specific signature ?

# GO Enrichment analysis 
# extract oligo gene blups
testGenes=oligo.stats$binarized.gene.blup[match(unique(oligo.stats[!oligo.stats$drop & !oligo.stats$dubious,]$GENEID), oligo.stats$GENEID)]
names(testGenes)=oligo.stats$GENEID[match(unique(oligo.stats[!oligo.stats$drop,]$GENEID), oligo.stats$GENEID)]
testGenes=sort(testGenes, decreasing=T)
GO.out=doGO(testGenes)
WriteXLS(GO.out, '/home/jbloom/Dropbox/Public/CoupledCRISPR/GO_gene_blups_ks_test.xls', SheetNames=names(GO.out))

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
plot(jitter(oligo.stats$hmm2[!oligo.stats$drop]), oligo.stats$binarized.oligo.blup[!oligo.stats$drop])
plot(jitter(oligo.stats$hmm2[!oligo.stats$drop]), oligo.stats$ALL.slope[!oligo.stats$drop])

osg=oligo.stats[!oligo.stats$drop , ]
osg=split(osg, osg$GENEID)
toleratePartialTrunc=sapply(osg, function(x) {
                nna=which(!is.na(x$hmm2) & x$ALL.cnt>5 )
                if(length(nna)==0) {FALSE} else {
                sum(x$hmm2[nna[1:2]]=='A')==2 & sum(x$binarized.oligo.blup[nna[1:2]]>0)==2  & sum(x$hmm2[nna[length(nna)]]=='D')==1 }
           })
names(which(toleratePartialTrunc))
allDead=sapply(osg, function(x) {
                nna=which(!is.na(x$hmm2))
                if(length(nna)==0) {FALSE }
                else {   y=x$hmm2[nna]
                      return(sum(y=='D')==length(y))
                    } 
              }
                )
names(which(allDead))
allAlive=sapply(osg, function(x) {
                nna=which(!is.na(x$hmm2))
                if(length(nna)==0) {FALSE }
                else {   y=x$hmm2[nna]
                      return(sum(y=='A')==length(y))
                    } 
              }
                )
names(which(allAlive))


GO.out.threshold2=doGO.threshold(names(which(allDead)), unique(oligo.stats$GENEID[!oligo.stats$drop]))


# some code to anotate the effects of stops in dubious orfs on the essential genes  ---------------
odub=oligo.stats[oligo.stats$dubious, ]
dubEffect = GRanges(seqnames=odub$chr, 
                        ranges=IRanges(start=odub$start, end=odub$end),
                        #strand=odub$codingStrand, 
                        strand=ifelse(odub$codingStrand=='+', '-', '+'),
                        id=odub$unique.Index,
                        gene=odub$GENEID
                        )
       #                 flag.syn=odub$flagsyn,

replacementSeq=DNAStringSet(ifelse(odub$codingStrand=='-',  'TCA', 'TGA'))
pdub=predictCoding(dubEffect,  txdb, sacCer3, replacementSeq, ignore.strand=FALSE)
pdub.df=data.frame(pdub)

binarized.mm.ng=glmer(slope.binarized~cnt+p_intercept+expt+binot+U6.Terminator+Score+(1|oligo), family=binomial(link="logit"), data=big.mm, verbose=T)
binarized.ranef.ng=ranef(binarized.mm.ng, cond=TRUE)
osng=binarized.ranef.ng$oligo[match(oligos$oligo, rownames(binarized.ranef.ng$oligo)),1]
pdub$blup=osng[pdub$id]
pdub$blupg=oligo.stats$binarized.oligo.blup[pdub$id]


boxplot(osgs.table[,1], osgs.table[,2])
a=cbind(osgs.table[,1], 1)
b=cbind(osgs.table[,2], 2)
cc=cbind(a,b)
apply(cc, 1, function(x) { segments(x[2],x[1], x[4], x[3]) })


osgs.table.bad=t(sapply(osgs, function(x) c(x$binarized.oligo.blup[x$dubious], x$binarized.oligo.blup[!x$dubious])))
t.test(osgs.table.bad[,1], osgs.table.bad[,2], paired=F)















gammALL.dd  =gamm4(slope.binarized~s(dist_from_CDS_end) +scale(cnt) +as.factor(expt)+ scale(p_intercept)+ scale(guide.GCcontent)+binot, random=~(1|oligo), family=binomial(link='logit'), data=big.mm[!big.mm$drop & big.mm$domain.downstream,]) 
gammALL.ndd =gamm4(slope.binarized~s(dist_from_CDS_end) +scale(cnt) +as.factor(expt)+ scale(p_intercept)+ scale(guide.GCcontent)+binot, random=~(1|oligo), family=binomial(link='logit'), data=big.mm[!big.mm$drop & !big.mm$domain.downstream,])

slopeAB=glmer(slope.binarized~
                cnt+
                as.factor(expt)+
                p_intercept+
                guide.GCcontent+
                dist_from
            
            _CDS_end+
                CDS_length+
                PAMvariantCNT+
                binot+
                lcd+
                evolvability+
                mean.aa.perfect.conserved+
                viable_annotation+
                end.conservation+
                domain.downstream+
                #spgf+
                (1|GENEID)+(1|oligo), data=b3, 
            family=binomial(link="logit"),
            control=glmerControl(optimizer="Nelder_Mead", optCtrl = list(maxfun = 1e6)) , verbose=T)
 pp=predict(slopeAB, type='response')
pp2=residuals(slopeAB) #, type='response')    
plot(b3$dist_from_CDS_end[match(names(pp), rownames(b3))], pp2, xlim=c(150,0))

plot(b3$dist_from_CDS_end, jitter(as.numeric(b3$hmm2)), xlim=c(150,0))

    #summary(slope5)
    gALL.anova=(Anova(slopeAB, type='III'))

    at=(anova(slopeAB, test='F'))
    length(predict(slopeA))
    at[order(at[,1], decreasing=T),]
ssc=summary(slopeAB)$coefficients

# prototype with slopes ... switch back to binarized if possible

s2=lmer(slope~  cnt+
                as.factor(expt)+
                p_intercept+
                guide.GCcontent+
                #dist_from_CDS_end+
                CDS_length+
                PAMvariantCNT+
                binot+
                #lcd+
                #evolvability+
                #mean.aa.perfect.conserved+
                #viable_annotation+
                #end.conservation+
                #domain.downstream+
                #spgf+
                (1|GENEID)+(1|oligo), data=b3, verbose=T)
s2a=(Anova(s2, type='III'))
s2resid=residuals(s2)
s2ranef=ranef(s2)

plot(jitter(b3$dist_from_CDS_end),b3$slope, xlim=c(150,0), col="#00000044",ylim=c(-.125, .05),
     xlab=c("AA distance from 3' end"), ylab='barcode slope')
plot(oligo.stats$dist_from_CDS_end[match(rownames(s2ranef$oligo), oligo.stats$oligo)], s2ranef$oligo[,1], xlim=c(1000,0), col="#00000055")
x2=cbind(oligo.stats$dist_from_CDS_end[match(rownames(s2ranef$oligo), oligo.stats$oligo)], s2ranef$oligo[,1])
colnames(x2)=c('distx', 'y')
x2=data.frame(x2)
l1=lm(y~distx, data=x2)
s3=segmented(l1, seg.Z=~distx)

boxplot(cut(x2$y~x2$distx, 100))

, col="#00000044",ylim=c(-.125, .05),
     xlab=c("AA distance from 3' end"), ylab='barcode slope')


plot(b3$dist_from_CDS_end[match(names(s2r), rownames(b3))],s2r, xlim=c(150,0))

sspline=smooth.spline(b3$dist_from_CDS_end[match(names(s2r), rownames(b3))],s2r)

# xx=cbind(x$dist_from_CDS_end[x$domain.downstream], x$binarized.oligo.blup.EssentialAll[x$domain.downstream])
    xx=na.omit(xx)
    data=xx
    sc=spline.cis(data,B=1000, alpha=.01)

#gAp=plot(gammALL$gam)
#gApD=plot(gammALL.dd$gam)
#gApnD=plot(gammALL.ndd$gam)

#dev.off()
#plot(gAp[[1]]$x, gAp[[1]]$fit, type='l', xlim=c(200,0))
#points(gAp[[1]]$x, gAp[[1]]$fit+gAp[[1]]$se, type='l', lty=2)
#points(gAp[[1]]$x, gAp[[1]]$fit-gAp[[1]]$se, type='l', lty=2)

# raw data and lines from regression model











###### NEW ###################################################3
#remove gene effect
#ibmm=glmer(slope.binarized~
#                cnt+
#                as.factor(expt)+
#                p_intercept+
#                dist_from_CDS_end+
#                CDS_length+
#                PAMvariantCNT+
#                binot+
#                (1|oligo),
#                 family=binomial(link="logit"), data=big.mm[!big.mm$drop,], verbose=T)
#rbmm=ranef(ibmm, cond=TRUE)
#om =oligo.stats

#dA=cbind(rbmm$oligo[,1],as.vector(attr(rbmm$oligo, 'postVar')) )
#rownames(dA)=rownames(rbmm$oligo)
#    colnames(dA)=c('b.oligo.blup', 'b.oligo.blup.postVar' )

# rewrite without using merge ()
#om= data.frame(om, dA[match(oligo.stats$oligo, rownames(dA)),], stringsAsFactors=F)

#og =split(om$b.oligo.blup, om$GENEID)
#ogs=sapply(og, function(x) sum(x<0, na.rm=T)/length(x))

#ogd=unique(om$GENEID[om$dubious])
#ogdd=unique(om$GENEID[!om$dubious])

#x11()
#par(mfrow=c(2,1))
#hist(ogs[ogd], breaks=100, main='dubious', xlim=c(0,1), xlab='fraction of oligos dead')
#hist(ogs[ogdd], breaks=100, main='non-dubious',xlim=c(0,1), xlab='fraction of oligos dead')

#x11()
#par(mfrow=c(2,1))
#hist(om$b.oligo.blup[!om$dubious], breaks=100, xlim=c(-3,4))
#hist(om$b.oligo.blup[om$dubious], breaks=100, xlim=c(-3,4))
#########################################################################

#average amount of conservation at specific positions 


