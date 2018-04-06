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
library(segmented)
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

out.base.dir='/media/jbloom/d1/coupled_CRISPR/Experiments/082316/'

# output here 
# dir.create(paste(out.base.dir, 'processed/RData/', sep=''))
# load oligos and annotation 
source('/media/jbloom/d1/coupled_CRISPR/code/accessory_functions.R')

# custom HMM code 
source('/media/jbloom/d1/coupled_CRISPR/code/mobsHMM.R')

seq.tables=readRDS(file = paste0(out.base.dir, 'processed/RData/seq.tables.RData'))
#giant.table=readRDS(file = paste0(out.base.dir, 'processed/RData/giant.table.RData'))
giant.table=readRDS(file = paste0(out.base.dir, 'processed/RData/giant.table.aligned2.RData'))

# minimum number of reads at t0
cutoff_t0=20
giant.table=addFilters(giant.table,cutoff_t0)

# Diagnostic Plots  --------------------------------------------------------------------------------------------------------------------
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
                #9:13,14, 47, 120)]
              x2= x[,c('WT_0', 'WT_24', 'WT_48', 'WT_72', 'WT_96', 'WT_72plate', 'matched.barcode')] 
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
   
    #low complexity domain indicator
    big.mm$lcd=as.numeric(ifelse(big.mm$low.complexity.downstream.cnt>0, 1, 0))
    # categroize as having off target gRNAs
    big.mm$binot=ifelse(big.mm$cnt.offtarget==1,0,1)
    
   # save(big.mm, file = paste0(out.base.dir, 'processed/RData/big.mm.RData'))
   # save(big.mm, file='/home/jbloom/Dropbox/Public/CoupledCRISPR/big.mm.RData')
    #load(paste0(out.base.dir, 'processed/RData/big.mm.RData'))
#--------------------------------------------------------------------------------- 
pdf('/home/jbloom/Dropbox/Public/CoupledCRISPR/SuppFig13_thetas.pdf', width=6, height=5)
hist(big.mm$slope, breaks=100, col='lightgrey', main='', xlab='barcoded plasmid theta')
abline(v=-0.025, col='red', lwd=2)
dev.off()

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

# Fit models with slopes --------------------------------------------------------------------------------------------------------------------------------

    # all data with gene term
    formula.input=as.formula(slope~cnt+p_intercept+expt+binot+U6.Terminator+Score+guide.GCcontent+PAMvariantCNT+(1|GENEID)+(1|oligo))
    oligo.stats=doLMER.oligo.gene(formula.input, data.in=big.mm, os, lab='')

    # no gene term
    formula.input=as.formula(slope~cnt+p_intercept+expt+binot+U6.Terminator+Score+guide.GCcontent+PAMvariantCNT+(1|oligo))
    oligo.stats=doLMER.oligo.gene(formula.input, data.in=big.mm, oligo.stats, lab='ng', doGene=FALSE)

    # essentials only
    formula.input=as.formula(slope~cnt+p_intercept+expt+binot+U6.Terminator+Score+guide.GCcontent+PAMvariantCNT+(1|GENEID)+(1|oligo))
    oligo.stats=doLMER.oligo.gene(formula.input, data.in=big.mm[!big.mm$drop,], oligo.stats, lab='EssentialAll')

    # WT only and no gene term
    formula.input=as.formula(slope~cnt+p_intercept+binot+U6.Terminator+Score+guide.GCcontent+PAMvariantCNT+(1|oligo))
    oligo.stats=doLMER.oligo.gene(formula.input, data.in=big.mm[ (big.mm$expt == 'WT'),], oligo.stats, lab='ng.WT', doGene=FALSE)

    # NMD only and no gene term
    formula.input=as.formula(slope~cnt+p_intercept+binot+U6.Terminator+Score+guide.GCcontent+PAMvariantCNT+(1|oligo))
    oligo.stats=doLMER.oligo.gene(formula.input, data.in=big.mm[ (big.mm$expt == 'NMD'),], oligo.stats, lab='ng.NMD', doGene=FALSE)

    # essentials WT only
    formula.input=as.formula(slope~cnt+p_intercept+binot+U6.Terminator+Score+guide.GCcontent+PAMvariantCNT+(1|GENEID)+(1|oligo))
    oligo.stats=doLMER.oligo.gene(formula.input, data.in=big.mm[!big.mm$drop & (big.mm$expt == 'WT'),], oligo.stats, lab='EssentialWT')

    # essentials NMD only
    formula.input=as.formula(slope~cnt+p_intercept+binot+U6.Terminator+Score+guide.GCcontent+PAMvariantCNT+(1|GENEID)+(1|oligo))
    oligo.stats=doLMER.oligo.gene(formula.input, data.in=big.mm[!big.mm$drop & (big.mm$expt == 'NMD'),], oligo.stats, lab='EssentialNMD')

    # essentials with downstream domain
    formula.input=as.formula(slope~cnt+p_intercept+binot+U6.Terminator+Score+guide.GCcontent+PAMvariantCNT+(1|GENEID)+(1|oligo))
    oligo.stats=doLMER.oligo.gene(formula.input, data.in=big.mm[!big.mm$drop & big.mm$domain.downstream,], oligo.stats, lab='Essential.DomainDownstream')

    # essentials with no downstream domain
    formula.input=as.formula(slope~cnt+p_intercept+binot+U6.Terminator+Score+guide.GCcontent+PAMvariantCNT+(1|GENEID)+(1|oligo))
    oligo.stats=doLMER.oligo.gene(formula.input, data.in=big.mm[!big.mm$drop & !big.mm$domain.downstream,], oligo.stats, lab='Essential.noDomainDownstream')

#-----------------------------------------------------------------------------------------------------------------------------------------------------------



# Fit models with binarized slopes -------------------------------------------------------------------------------------------------------------------------

    # all data with gene term
    formula.input=as.formula(slope.binarized~cnt+p_intercept+expt+binot+U6.Terminator+Score+guide.GCcontent+PAMvariantCNT+(1|GENEID)+(1|oligo))
    oligo.stats=doGLMER.oligo.gene(formula.input, data.in=big.mm, oligo.stats, lab='')

    # no gene term
    formula.input=as.formula(slope.binarized~cnt+p_intercept+expt+binot+U6.Terminator+Score+guide.GCcontent+PAMvariantCNT+(1|oligo))
    oligo.stats=doGLMER.oligo.gene(formula.input, data.in=big.mm, oligo.stats, lab='ng', doGene=FALSE)
    
    # essentials only
    formula.input=as.formula(slope.binarized~cnt+p_intercept+expt+binot+U6.Terminator+Score+guide.GCcontent+PAMvariantCNT+(1|GENEID)+(1|oligo))
    oligo.stats=doGLMER.oligo.gene(formula.input, data.in=big.mm[!big.mm$drop,], oligo.stats, lab='EssentialAll')
    
    # WT only and no gene term
    formula.input=as.formula(slope.binarized~cnt+p_intercept+binot+U6.Terminator+Score+guide.GCcontent+PAMvariantCNT+(1|oligo))
    oligo.stats=doGLMER.oligo.gene(formula.input, data.in=big.mm[ (big.mm$expt == 'WT'),], oligo.stats, lab='ng.WT', doGene=FALSE)

    # NMD only and no gene term
    formula.input=as.formula(slope.binarized~cnt+p_intercept+binot+U6.Terminator+Score+guide.GCcontent+PAMvariantCNT+(1|oligo))
    oligo.stats=doGLMER.oligo.gene(formula.input, data.in=big.mm[ (big.mm$expt == 'NMD'),], oligo.stats, lab='ng.NMD', doGene=FALSE)

    # essentials WT only
    formula.input=as.formula(slope.binarized~cnt+p_intercept+binot+U6.Terminator+Score+guide.GCcontent+PAMvariantCNT+(1|GENEID)+(1|oligo))
    oligo.stats=doGLMER.oligo.gene(formula.input, data.in=big.mm[!big.mm$drop & (big.mm$expt == 'WT'),], oligo.stats, lab='EssentialWT')

    # essentials NMD only
    formula.input=as.formula(slope.binarized~cnt+p_intercept+binot+U6.Terminator+Score+guide.GCcontent+PAMvariantCNT+(1|GENEID)+(1|oligo))
    oligo.stats=doGLMER.oligo.gene(formula.input, data.in=big.mm[!big.mm$drop & (big.mm$expt == 'NMD'),], oligo.stats, lab='EssentialNMD')

    # essentials with downstream domain
    formula.input=as.formula(slope.binarized~cnt+p_intercept+binot+U6.Terminator+Score+guide.GCcontent+PAMvariantCNT+(1|GENEID)+(1|oligo))
    oligo.stats=doGLMER.oligo.gene(formula.input, data.in=big.mm[!big.mm$drop & big.mm$domain.downstream,], oligo.stats, lab='Essential.DomainDownstream')

    # essentials with no downstream domain
    formula.input=as.formula(slope.binarized~cnt+p_intercept+binot+U6.Terminator+Score+guide.GCcontent+PAMvariantCNT+(1|GENEID)+(1|oligo))
    oligo.stats=doGLMER.oligo.gene(formula.input, data.in=big.mm[!big.mm$drop & !big.mm$domain.downstream,], oligo.stats, lab='Essential.noDomainDownstream')


    oligo.stats=oligo.stats[order(oligo.stats$unique.Index),]
#-----------------------------------------------------------------------------------------------------------------------------------------------------------


# === HMM ================================================================================================
    # do HMM 
    hmm.out=doHMM(big.mm, oligo.stats)
    oligo.stats = data.frame(oligo.stats, hmm.out[match(oligos$oligo, hmm.out$name),])
 
    acnt=sapply(split(oligo.stats$hmm2, oligo.stats$GENEID), function(x) sum(x=='A',na.rm=T))
    afrac=sapply(split(oligo.stats$hmm2, oligo.stats$GENEID), function(x) sum(x=='A',na.rm=T)/sum(!is.na(x)))
    atable=cbind(acnt, afrac)
    colnames(atable)=c('hmm.count.alive', 'hmm.frac.alive')
    rm(acnt, afrac)
    oligo.stats= data.frame(oligo.stats, atable[match(oligo.stats$GENEID, rownames(atable)),], stringsAsFactors=F)
    oligo.stats=oligo.stats[order(oligo.stats$unique.Index),]

    # SEND  supp fig 6
    pdf(file='/home/jbloom/Dropbox/Public/CoupledCRISPR/Figures/SuppFig6_HMM_fraction_per_gene.pdf', width=5, height=5)
    hist(acnt[!is.na(afrac)], breaks=50, xlab='number of PTCs tolerated', ylab='genes', main='', xaxt='n')
    axis(1, c(0:10))
    dev.off()
    #-------------------------------------------------------------------
    plot.dir='/media/jbloom/d1/coupled_CRISPR/plots/HMM_v10/'
    #makeHMMplots(plot.dir, oligo.stats, big.mm, conservation, dgsplit)
    
   # save(oligo.stats, file=paste0(out.base.dir, 'processed/RData/oligo.stats.EssentialBlups.RData'))
   # save(oligo.stats, file='/home/jbloom/Dropbox/Public/CoupledCRISPR/oligo.stats.EssentialBlups.RData')
    WriteXLS(oligo.stats, '/home/jbloom/Dropbox/Public/CoupledCRISPR/oligo_stats.xls', row.names=F)

    #load(paste0(out.base.dir, 'processed/RData/oligo.stats.EssentialBlups.RData'))
# ====================================================================================================



#============================================================================================

splic.genes=unlist(strsplit(GO.out[[1]][1,]$genes.in.set.sorted, ' '))
cat.genes=unlist(strsplit(GO.out[[2]][1,]$genes.in.set.sorted, ' '))
ostatG=oligo.stats[match(unique(oligo.stats$GENEID), oligo.stats$GENEID),]
dff= ostatG[,c("GENEID","binarized.gene.blup","binarized.gene.blup.EssentialAll", "hmm.count.alive","hmm.frac.alive")]
WriteXLS(dff,'/home/jbloom/Dropbox/Public/CoupledCRISPR/gene_stats.xls', row.names=F)


spgf=ostatGoligo.stats$GENEID %in% splic.genes
cgf=ostatG$GENEID %in% cat.genes
# SEND putative supplementary figure for GO enrichment
pdf(file='/home/jbloom/Dropbox/Public/CoupledCRISPR/Figures/SuppFig_GO_enrichment.pdf', width=12, height=10)
par(mfrow=c(1,2))
stripchart(ostatG$binarized.gene.blup.EssentialAll~spgf, vertical=T, method='jitter', pch=20, col='#00000055', 
           main='RNA splicing, via transesterification reactions GO:0000375', xlim=c(.5,2.5), ylab='gene PTC tolerance score', group.names=rev(c('genes in category', 'other genes')) )
boxplot(ostatG$binarized.gene.blup.EssentialAll~spgf, add=T,outline=F, xaxt='n')
stripchart(ostatG$binarized.gene.blup.EssentialAll~cgf, vertical=T, method='jitter', pch=20, col='#00000055', 
           main='catalytic activity   GO:0003824', xlim=c(.5,2.5), ylab='gene PTC tolerance score', group.names=rev(c('genes in category', 'other genes')))
boxplot(ostatG$binarized.gene.blup.EssentialAll~cgf,add=T ,outline=F, xaxt='n')
dev.off()    
    
    
    # BIG model with all----------------------------------------------------------------------------------- 
b3=big.mm[!big.mm$drop , ] # & (big.mm$expt == 'WT') , ]#& big.mm$dist_from_CDS_end<300,]
b3$expt=as.factor(b3$expt)
b3$oligo=as.factor(b3$oligo)
b3$domain.downstream=as.factor(b3$domain.downstream)
# scale(H0_w)+
splic.genes=unlist(strsplit(GO.out[[1]][1,]$genes.in.set.sorted, ' '))
spgf=b3$GENEID %in% splic.genes

cat.genes=unlist(strsplit(GO.out[[2]][1,]$genes.in.set.sorted, ' '))
cgf=b3$GENEID %in% cat.genes
fractionOfGene=b3$dist_from_CDS_end/b3$CDS_length

s2=glmer(slope.binarized~cnt+expt+p_intercept+guide.GCcontent+dist_from_CDS_end+CDS_length+PAMvariantCNT+
                binot+U6.Terminator+Score+lcd+evolvability+mean.aa.perfect.conserved+viable_annotation+
                end.conservation+domain.downstream+fractionOfGene+(1|GENEID)+(1|oligo), data=b3, verbose=T,  family=binomial(logit), 
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

#goID=gt[3,'GO.ID']
#plot(showGroupDensity(GOData, goID, ranks=T) )



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
#on2=oligo.stats[match(unique(oligo.stats$GENEID), oligo.stats$GENEID),]

# change axis labels (PTC tolerance to PTC tolerance score)
# change (codons from C' terminal) to (PTC distance from 3' end (codons) )


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
    
    Ws=split(oligo.stats$binarized.oligo.blup.ng.WT, oligo.stats$dubious)
    names(Ws)=c('WT Essential', 'WT Dubious')
    Ns=split(oligo.stats$binarized.oligo.blup.ng.NMD, oligo.stats$dubious)
    names(Ns)=c('NMD Essential', 'NMD Dubious')

    Blist=(c(Ws,Ns))
    Blist=Blist[c(1,3,2,4)]

    par(mfrow=c(2,1))
    stripchart(Blist, vertical=T, method='jitter', col='#00000012', pch=19, at=c(1,1.5,2.25,2.75), xlim=c(.75,3.25),ylab='PTC tolerance')
    boxplot(Blist, add=T, xlab='', names=rep('',4), ylab='', outline=F, xaxt='n',
            at=c(1.25,1.75, 2.5, 3), 
            boxwex=.15) 

    Blist2=split(oligo.stats$binarized.oligo.blup.ng,oligo.stats$dubious)
    stripchart(Blist2, vertical=T, method='jitter', col='#00000012', pch=19, at=c(1,2), xlim=c(.75,2.5),ylab='PTC tolerance')
    boxplot(Blist2, add=T, xlab='', names=rep('',2), ylab='', outline=F, xaxt='n',
            at=c(1.25, 2.25), 
            boxwex=.25) 



    test= stripchart(oligo.stats$binarized.oligo.blup.ng~oligo.stats$dubious, vertical=T, method='jitter', col='#00000012', pch=19,
               group.names=c('essential', 'dubious'),
               ylab='PTC tolerance', xlim=c(0.5,2.5) )
    boxplot(oligo.stats$binarized.oligo.blup.ng~oligo.stats$dubious, add=T, xlab='', names='', ylab='', outline=F) 
    dev.off()

    dup.guides=oligo.stats$guide[which(duplicated(oligo.stats$guide))]
    idg= oligo.stats$guide %in% dup.guides
     
    osgs=split(oligo.stats[idg,], oligo.stats[idg,]$guide)
    osgs.table=t(sapply(osgs, function(x) c(x$binarized.oligo.blup.ng[x$dubious], x$binarized.oligo.blup.ng[!x$dubious])))
    osgs.table.id=t(sapply(osgs, function(x) c(x$unique.Index[x$dubious], x$unique.Index[!x$dubious])))
    colnames(osgs.table)=c('dubious', 'essential')
    colnames(osgs.table.id)=c('dubious', 'essential')

   
    odub=oligo.stats[oligo.stats$dubious, ] #[as.numeric(osgs.table.id[,'dubious']),]
    dubEffect = GRanges(seqnames=odub$chr, 
                        ranges=IRanges(start=odub$start, end=odub$end),
                        #strand=odub$codingStrand, 
                        strand=ifelse(odub$codingStrand=='+', '-', '+'),
                        unique.Index=odub$unique.Index,
                        id=rownames(odub), 
                        gene=odub$GENEID
                        )

    replacementSeq=DNAStringSet(ifelse(odub$codingStrand=='-',  'TCA', 'TGA'))
    pdub=predictCoding(dubEffect,  txdb, sacCer3, replacementSeq, ignore.strand=FALSE)
   
    prov.effects=list()
    for(i in 1:length(pdub)) {
    goi=pdub$GENEID[i]
    loi=pdub$PROTEINLOC[[i]]
    varoi=s2c(as.character(pdub$VARAA[i]))
    idoi=pdub$unique.Index[i]
    refoi=s2c(as.character(pdub$REFAA[i]))
    prot.seq.rownames=rownames(provean_scores[[goi]])

    if(paste(prot.seq.rownames[loi], collapse='')==paste(refoi, collapse='')) {
        peff=try({provean_scores[[goi]][cbind(loi,match(varoi, colnames(provean_scores[[goi]])  ))]})
        if(class(peff)=='try-error') { prov.effects[[idoi]]= NA } 
        else {    prov.effects[[idoi]]=provean_scores[[goi]][cbind(loi,match(varoi, colnames(provean_scores[[goi]])))] }
    } else {prov.effects[[idoi]]=NA }
    }

    max.prov.effect=sapply(prov.effects, function(x) x[which.min((x))] )
    pscores=sapply(max.prov.effect, function(x) if(length(x)==0) {return(NA)} else {x } )

    plot(oligo.stats$binarized.oligo.blup.ng[1:length(pscores)], pscores, ylab='provean score', xlab='PTC tolerance')
    cor.test(oligo.stats$binarized.oligo.blup.ng[1:length(pscores)], pscores, method='spearman')

    plot(oligo.stats$slope.oligo.blup.ng[1:length(pscores)], pscores)
    cor.test(oligo.stats$slope.oligo.blup.ng[1:length(pscores)], pscores, method='spearman')
     



      grantham=read.delim('/media/jbloom/d1/coupled_CRISPR/Reference/grantham.matrix.csv', header=T, sep='\t')
    rownames(grantham)=a(rownames(grantham))
    colnames(grantham)=a(colnames(grantham))
    
    granth.effects=list() 
    for(i in 1:length(pdub)) {
    goi=pdub$GENEID[i]
    loi=pdub$PROTEINLOC[[i]]
    varoi=s2c(as.character(pdub$VARAA[i]))
    refoi=s2c(as.character(pdub$REFAA[i]))
    idoi=pdub$unique.Index[i]
 
    ec=grep('\\*', refoi)
    if(length(ec)>0 ) {
        varoi=varoi[-ec]
        refoi=refoi[-ec]
    }
    vrlookup=cbind(refoi, varoi)
    granth.effects[[idoi]]=grantham[vrlookup]
    }
    max.granth.effect=sapply(granth.effects, function(x) x[which.max((x))] )
    gscores=sapply(max.granth.effect, function(x) if(length(x)==0) {return(NA)} else {x } )
    gscores[gscores==0]=NA
    plot(oligo.stats$binarized.oligo.blup.ng[1:length(gscores)], gscores, ylab='grantham score', xlab='PTC tolerance')
    cor.test(oligo.stats$binarized.oligo.blup.ng[1:length(gscores)], gscores, method='spearman')






    oss=oligo.stats[1:length(pscores),]
    oss$pscores=pscores

    formula.input=as.formula(slope~cnt+p_intercept+expt+binot+U6.Terminator+Score+guide.GCcontent+PAMvariantCNT+(1|GENEID)+(1|oligo))
   
    b2m=big.mm
    b2m$pscores=NA
    b2m$pscores=pscores[!is.na(pscores)][match(b2m$unique.Index, which(!is.na(pscores)))]
    b2l=(glmer(slope.binarized~cnt+p_intercept+expt+binot+U6.Terminator+Score+guide.GCcontent+PAMvariantCNT+dist_from_CDS_end+(1|oligo), data=b2m, family=binomial('logit')))
    Anova(b2l, type='III')
    
    
    
    drop1(lm(slope.oligo.blup.ng~U6.Terminator+Score+guide.GCcontent+PAMvariantCNT+dist_from_CDS_end+pscores, data=oss), test='Chisq')
   
   
    #, method='spearman')
   
    pscores.sub=pscores[osgs.table.id[,'dubious']]

    color.gradient <- function(x, colors=c("red", "blue"), colsteps=100) {
        return( colorRampPalette(colors, bias=.5) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
    }
    pcolors.sub=rep(NA, length(pscores.sub))
    pcolors.sub[!is.na(pscores.sub)]=color.gradient(pscores.sub[!is.na(pscores.sub)])

    plot(rep(1, nrow(osgs.table)), osgs.table[,2], xlim=c(-0,3), ylim=c(-2,4.2),col =pcolors)
    points(rep(2, nrow(osgs.table)), osgs.table[,1], xlim=c(-0,3), ylim=c(-2, 4.2), col =pcolors)
    segments(rep(1, nrow(osgs.table)), osgs.table[,2], rep(2, nrow(osgs.table)),  osgs.table[,1], col =pcolors)
        

    pdf(file='/home/jbloom/Dropbox/Public/CoupledCRISPR/Figures/sFigure2_matched_gRNA.pdf', width=8, height=8)
    plot(osgs.table, xlim=c(-2,4), ylim=c(-2,4), ylab='essential PTC tolerance', xlab='dubious PTC tolerance' )
    abline(0,1)
    dev.off()
    t.test(osgs.table[,1], osgs.table[,2], paired=T)
    load('/data/Databases/provean_scores.RData')

#------------------------------------------------------------------------------------------------------------




# Figure 3 -----------------------------------------------------------------------------------------------
endcons.thresh=.8306452
column="binarized.oligo.blup.EssentialAll"
#column="binarized.oligo.blup.Essential.DomainDownstream"
#column="binarized.oligo.blup.Essential.noDomainDownstream"
column="binarized.oligo.blup.EssentialWT"
column="binarized.oligo.blup.EssentialNMD"
pdf(file='/home/jbloom/Dropbox/Public/CoupledCRISPR/Figures/Figure3_150.pdf', width=7.5, height=5)
pdf(file='/home/jbloom/Dropbox/Public/CoupledCRISPR/Figures/Figure3_150_DomaninvNoDomain.pdf', width=7.5, height=10)
pdf(file='/home/jbloom/Dropbox/Public/CoupledCRISPR/Figures/Figure3_150_EndConsHigh_v_EndConsLow.pdf', width=7.5, height=10)

pdf(file='/home/jbloom/Dropbox/Public/CoupledCRISPR/Figures/Figure3_150_WTtopvNMDbot.pdf', width=8, height=10)
par(mfrow=c(2,1))
column="binarized.oligo.blup.EssentialWT"
plot(oligo.stats$dist_from_CDS_end, oligo.stats[,column], xlim=c(150,0),xaxt='n', ylab='PTC tolerance', xlab='AA distance from C-terminal', col='#00000022',cex=.75 ,pch=20 ,main='WT', ylim=c(-3,3.5))
axis(1, at=rev(seq(0,150,5)))
x2=na.omit(data.frame(distx=oligo.stats$dist_from_CDS_end, y=oligo.stats[,column]))
x3=x2[x2$distx<500,]
s3=segmented(lm(y~distx, data=x3), seg.Z=~distx, control=seg.control(n.boot=100))
seg.fit=cbind(x3$distx, predict(s3))
seg.fit=seg.fit[order(seg.fit[,1]),]
points(seg.fit[,1], seg.fit[,2], col='blue', type='l', lwd=4)
segments(confint.segmented(s3)$distx[c(2,3)], c(-.5,-.5) , confint.segmented(s3)$distx[c(2,3)] , c(.5,.5), col='blue', lwd=2)
column="binarized.oligo.blup.EssentialNMD"
plot(oligo.stats$dist_from_CDS_end, oligo.stats[,column], xlim=c(150,0),xaxt='n', ylab='PTC tolerance', xlab='AA distance from C-terminal', col='#00000022',cex=.75 ,pch=20 ,main='nmd(del)', ylim=c(-3,3.5))
axis(1, at=rev(seq(0,150,5)))
x2=na.omit(data.frame(distx=oligo.stats$dist_from_CDS_end, y=oligo.stats[,column]))
x3=x2[x2$distx<500,]
s3=segmented(lm(y~distx, data=x3), seg.Z=~distx, control=seg.control(n.boot=100))
seg.fit=cbind(x3$distx, predict(s3))
seg.fit=seg.fit[order(seg.fit[,1]),]
points(seg.fit[,1], seg.fit[,2], col='blue', type='l', lwd=4)
segments(confint.segmented(s3)$distx[c(2,3)], c(-.5,-.5) , confint.segmented(s3)$distx[c(2,3)] , c(.5,.5), col='blue', lwd=2)
dev.off()

#(A) PTCs affecting less conserved sequence and disrupting an annotated protein domain.
#(B) PTCs affecting more conserved sequence and disrupting an annotated protein domain. 
#(c) PTCs affecting less conserved sequence and not disrupting an annotated protein domain. 
#(D) PTCs affecting more conserved sequence and not disrupting an annotated protein domain.
pdf(file='/home/jbloom/Dropbox/Public/CoupledCRISPR/Figures/SuppFigureDistance.pdf', width=15, height=10)
par(mfrow=c(2,2))
column="binarized.oligo.blup.EssentialAll"
plot(oligo.stats$dist_from_CDS_end[oligo.stats$end.conservation<=endcons.thresh & oligo.stats$domain.downstream], oligo.stats[,column][oligo.stats$end.conservation<=endcons.thresh & oligo.stats$domain.downstream], 
                                                                xlim=c(150,0),xaxt='n', ylab='PTC tolerance', xlab='AA distance from C-terminal', col='#00000044',cex=.75 ,pch=20 , ylim=c(-3,4), main='less conserved, disrupting a domain')
axis(1, at=rev(seq(0,150,5)))

x2=na.omit(data.frame(distx=oligo.stats$dist_from_CDS_end[oligo.stats$end.conservation<=endcons.thresh & oligo.stats$domain.downstream], y=oligo.stats[,column][oligo.stats$end.conservation<=endcons.thresh & oligo.stats$domain.downstream]))
x3=x2[x2$distx<500,]
s3=segmented(lm(y~distx, data=x3), seg.Z=~distx, control=seg.control(n.boot=100))
seg.fit=cbind(x3$distx, predict(s3))
seg.fit=seg.fit[order(seg.fit[,1]),]
points(seg.fit[,1], seg.fit[,2], col='blue', type='l', lwd=4)
#segments(confint.segmented(s3)$distx[c(2,3)], c(-.5,-.5) , confint.segmented(s3)$distx[c(2,3)] , c(.5,.5), col='blue', lwd=2)


plot(oligo.stats$dist_from_CDS_end[oligo.stats$end.conservation>endcons.thresh & oligo.stats$domain.downstream], oligo.stats[,column][oligo.stats$end.conservation>endcons.thresh & oligo.stats$domain.downstream], 
                                                                xlim=c(150,0),xaxt='n', ylab='PTC tolerance', xlab='AA distance from C-terminal', col='#00000044',cex=.75 ,pch=20 , ylim=c(-3,4),main='more conserved, disrupting a domain')
axis(1, at=rev(seq(0,150,5)))

x2=na.omit(data.frame(distx=oligo.stats$dist_from_CDS_end[oligo.stats$end.conservation>endcons.thresh & oligo.stats$domain.downstream], y=oligo.stats[,column][oligo.stats$end.conservation>endcons.thresh & oligo.stats$domain.downstream]))
x3=x2[x2$distx<500,]
s3=segmented(lm(y~distx, data=x3), seg.Z=~distx, control=seg.control(n.boot=100))
seg.fit=cbind(x3$distx, predict(s3))
seg.fit=seg.fit[order(seg.fit[,1]),]
points(seg.fit[,1], seg.fit[,2], col='blue', type='l', lwd=4)
#segments(confint.segmented(s3)$distx[c(2,3)], c(-.5,-.5) , confint.segmented(s3)$distx[c(2,3)] , c(.5,.5), col='blue', lwd=2)

plot(oligo.stats$dist_from_CDS_end[oligo.stats$end.conservation<=endcons.thresh & !oligo.stats$domain.downstream], oligo.stats[,column][oligo.stats$end.conservation<=endcons.thresh & !oligo.stats$domain.downstream], 
                                                                xlim=c(150,0),xaxt='n', ylab='PTC tolerance', xlab='AA distance from C-terminal', col='#00000044',cex=.75 ,pch=20, ylim=c(-3,4) ,main='less conserved, not disrupting a domain')
axis(1, at=rev(seq(0,150,5)))

x2=na.omit(data.frame(distx=oligo.stats$dist_from_CDS_end[oligo.stats$end.conservation<=endcons.thresh & !oligo.stats$domain.downstream], y=oligo.stats[,column][oligo.stats$end.conservation<=endcons.thresh & !oligo.stats$domain.downstream]))
x3=x2[x2$distx<500,]
s3=segmented(lm(y~distx, data=x3), seg.Z=~distx, control=seg.control(n.boot=100))
seg.fit=cbind(x3$distx, predict(s3))
seg.fit=seg.fit[order(seg.fit[,1]),]
points(seg.fit[,1], seg.fit[,2], col='blue', type='l', lwd=4)
#segments(confint.segmented(s3)$distx[c(2,3)], c(-.5,-.5) , confint.segmented(s3)$distx[c(2,3)] , c(.5,.5), col='blue', lwd=2)

plot(oligo.stats$dist_from_CDS_end[oligo.stats$end.conservation>endcons.thresh & !oligo.stats$domain.downstream], oligo.stats[,column][oligo.stats$end.conservation>endcons.thresh & !oligo.stats$domain.downstream], 
                                                                xlim=c(150,0),xaxt='n', ylab='PTC tolerance', xlab='AA distance from C-terminal', col='#00000044',cex=.75 ,pch=20, ylim=c(-3,4) ,main='more conserved, not disrupting a domain')
axis(1, at=rev(seq(0,150,5)))

x2=na.omit(data.frame(distx=oligo.stats$dist_from_CDS_end[oligo.stats$end.conservation>endcons.thresh & !oligo.stats$domain.downstream], y=oligo.stats[,column][oligo.stats$end.conservation>endcons.thresh & !oligo.stats$domain.downstream]))
x3=x2[x2$distx<500,]
s3=segmented(lm(y~distx, data=x3), seg.Z=~distx, control=seg.control(n.boot=100))
seg.fit=cbind(x3$distx, predict(s3))
seg.fit=seg.fit[order(seg.fit[,1]),]
points(seg.fit[,1], seg.fit[,2], col='blue', type='l', lwd=4)
#segments(confint.segmented(s3)$distx[c(2,3)], c(-.5,-.5) , confint.segmented(s3)$distx[c(2,3)] , c(.5,.5), col='blue', lwd=2)
dev.off()

# by gene length
clc=cut(oligo.stats$CDS_length, c(0, 200, 400, 6000))
column="binarized.oligo.blup.EssentialAll"
par(mfrow=c(3,1))
for(i in 1:3) {
plot(oligo.stats$dist_from_CDS_end[clc==levels(clc)[i]], oligo.stats[,column][clc==levels(clc)[i]], 
                                                                xlim=c(150,0),xaxt='n', ylab='PTC tolerance', xlab='AA distance from C-terminal', col='#00000044',cex=.75 ,pch=20 , ylim=c(-3,4), main=levels(clc)[i])
axis(1, at=rev(seq(0,150,5)))
x2=na.omit(data.frame(distx=oligo.stats$dist_from_CDS_end[clc==levels(clc)[i]], y=oligo.stats[,column][clc==levels(clc)[i]]))
x3=x2[x2$distx<500,]
s3=segmented(lm(y~distx, data=x3), seg.Z=~distx, control=seg.control(n.boot=100))
seg.fit=cbind(x3$distx, predict(s3))
seg.fit=seg.fit[order(seg.fit[,1]),]
points(seg.fit[,1], seg.fit[,2], col='blue', type='l', lwd=4)
segments(confint.segmented(s3)$distx[c(2,3)], c(-.5,-.5) , confint.segmented(s3)$distx[c(2,3)] , c(.5,.5), col='blue', lwd=2)

}



plot(oligo.stats$dist_from_CDS_end[oligo.stats$domain.downstream], oligo.stats[,column][oligo.stats$domain.downstream], xlim=c(150,0),xaxt='n', ylab='PTC tolerance', xlab='AA distance from C-terminal', col='#00000022',cex=.75 ,pch=20 ,main='')
plot(oligo.stats$dist_from_CDS_end[!oligo.stats$domain.downstream], oligo.stats[,column][!oligo.stats$domain.downstream], xlim=c(150,0),xaxt='n', ylab='PTC tolerance', xlab='AA distance from C-terminal', col='#00000022',cex=.75 ,pch=20 ,main='')
plot(oligo.stats$dist_from_CDS_end[oligo.stats$end.conservation>endcons.thresh], oligo.stats[,column][oligo.stats$end.conservation>endcons.thresh], xlim=c(150,0),xaxt='n', ylab='PTC tolerance', xlab='AA distance from C-terminal', col='#00000022',cex=.75 ,pch=20 ,main='')
plot(oligo.stats$dist_from_CDS_end[oligo.stats$end.conservation<=endcons.thresh], oligo.stats[,column][oligo.stats$end.conservation<=endcons.thresh], xlim=c(150,0),xaxt='n', ylab='PTC tolerance', xlab='AA distance from C-terminal', col='#00000022',cex=.75 ,pch=20 ,main='')
x2=na.omit(data.frame(distx=oligo.stats$dist_from_CDS_end, y=oligo.stats[,column]))
x2=na.omit(data.frame(distx=oligo.stats$dist_from_CDS_end[oligo.stats$domain.downstream], y=oligo.stats[,column][oligo.stats$domain.downstream]))
x2=na.omit(data.frame(distx=oligo.stats$dist_from_CDS_end[!oligo.stats$domain.downstream], y=oligo.stats[,column][!oligo.stats$domain.downstream]))
x2=na.omit(data.frame(distx=oligo.stats$dist_from_CDS_end[oligo.stats$end.conservation>endcons.thresh], y=oligo.stats[,column][oligo.stats$end.conservation>endcons.thresh]))
x2=na.omit(data.frame(distx=oligo.stats$dist_from_CDS_end[oligo.stats$end.conservation<=endcons.thresh], y=oligo.stats[,column][oligo.stats$end.conservation<=endcons.thresh]))
x3=x2[x2$distx<500,]
s3=segmented(lm(y~distx, data=x3), seg.Z=~distx, control=seg.control(n.boot=100))
seg.fit=cbind(x3$distx, predict(s3))
seg.fit=seg.fit[order(seg.fit[,1]),]
points(seg.fit[,1], seg.fit[,2], col='blue', type='l', lwd=4)
segments(confint.segmented(s3)$distx[c(2,3)], c(-.5,-.5) , confint.segmented(s3)$distx[c(2,3)] , c(.5,.5), col='blue', lwd=2)
         
#
ods=split(oligo.stats[,column], oligo.stats$dist_from_CDS_end)
odsq=t(sapply(ods, function(x) quantile(x, c(.25,.5,.75), na.rm=T)))
no=as.numeric(rownames(odsq))
segments(no-.5, odsq[,1], no+.5, odsq[,1], lwd=2)
segments(no-.5, odsq[,2], no+.5, odsq[,2], col='red', lwd=2)
segments(no-.5, odsq[,3], no+.5, odsq[,3], lwd=2)
     
dev.off()
#---------------------------------------------------------------------------------------------------------


save.image('/media/jbloom/d1/coupled_CRISPR/cc052417.RData')













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

pdf(file='/home/jbloom/Dropbox/Public/CoupledCRISPR/Figures/suppHMM_PTC_level.pdf', width=10, height=5, useDingbats=F)
plot(as.numeric(names(se2.frac)), se2.frac, xlim=c(150,0), cex=2, pch=20, ylab='fraction of PTCs tolerated', xlab='AA distance from C-terminal')
# cex=sqrt(sapply(se2, length)/10))

smb=matrix(0, length(se2.frac), 100)
for(i in 1:100) {
smb[,i]=sapply(se2, function(x) {y=sample(x, replace=T);  return(sum(y)/length(y)) })
}
segments(as.numeric(names(se2.frac)),
         apply(smb,1, quantile, .05),
         as.numeric(names(se2.frac)),
         apply(smb,1, quantile, .95), col='black')
dev.off()

##### data subset without domain downstream
par(mfrow=c(2,1))
e3=cbind(oligo.stats$hmm2, oligo.stats$dist_from_CDS_end)[!oligo.stats$domain.downstream,]
e3[,1]=as.numeric(e3[,1])-1
e3=na.omit(e3)
se3=split(e3[,1], e3[,2])
se3.frac=sapply(se3, function(x) sum(x)/length(x))
plot(as.numeric(names(se3.frac)), se3.frac, xlim=c(150,0), cex=sqrt(sapply(se3, length)/20), 
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

plot(as.numeric(names(se4.frac)), se4.frac, xlim=c(150,0), cex=sqrt(sapply(se4, length)/20), 
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
        
        

