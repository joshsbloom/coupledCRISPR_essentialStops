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
#library(msm)
#library(afex)
#library(car)

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

slope.mm=lmer(slope~cnt+p_intercept+expt+binot+U6.Terminator+Score+(1|GENEID)+(1|oligo), data=big.mm, REML=TRUE)
slope.ranef=ranef(slope.mm, cond=TRUE)

gb=cbind(slope.ranef$GENEID[,1], as.vector(attr(slope.ranef$GENEID, 'postVar')))
colnames(gb)=c('slope.gene.blup', 'slope.gene.blup.postVar')
rownames(gb)=rownames(slope.ranef$GENEID)

# corrected call
#binarized.feff=glmer(slope.binarized~cnt+p_intercept+expt+binot+U6.Terminator+Score, family=binomial(link="logit"), data=big.mm, verbose=T)
#reclassify
#mod=glm(slope.binarized~cnt+p_intercept+binot+U6.Terminator+Score, family=binomial(link="logit"), data=big.mm)
#mod=lm(slope~cnt+p_intercept+binot+U6.Terminator+Score, data=big.mm)
#big.mm$slope.corrected=residuals(mod)
#plot(big.mm$slope, residuals(mod), col='#00000010')
#big.mm$slope.binarized.corrected=ifelse(residuals(mod)>(0.032), 1,0)
#test1=glmer(slope.binarized.corrected~cnt+p_intercept+expt+binot+U6.Terminator+Score+dubious+(1|GENEID)+(1|oligo), family=binomial(link="logit"), data=big.mm, verbose=T)
#test2=glmer(slope.binarized~cnt+p_intercept+expt+binot+U6.Terminator+Score+dubious+(1|GENEID)+(1|oligo), family=binomial(link="logit"), data=big.mm, verbose=T)

# modified !!!!!!!
binarized.mm=glmer(slope.binarized~cnt+p_intercept+expt+binot+U6.Terminator+Score+(1|GENEID)+(1|oligo), family=binomial(link="logit"), data=big.mm, verbose=T)
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

#save.image(paste0(out.base.dir, 'processed/RData/workspace03092017.RData'))

#  HMM ================================================================================================
    # do HMM 
    hmm.out=doHMM(big.mm, oligo.stats)
    oligo.stats = data.frame(oligo.stats, hmm.out[match(oligos$oligo, hmm.out$name),])
    #hmm.out=doHMM(big.mm, oligo.stats)
    #save(hmm.out, file=paste0(out.base.dir, 'processed/RData/hmm.out.030917.RData'))
    #load(paste0(out.base.dir, 'processed/RData/hmm.out.032317.RData'))

    # integrate data from hmm 2 state into oligo table 
    #so=sapply(hmm.out, function(x) x$states.out)
    #
    #dso=data.frame(GENEID=rep(names(so), sapply(so, length )), 
    #               dist_from_CDS_end=as.vector(unlist(sapply(so, names))),
    #               hmm2=as.vector(unlist(so)))
    #dso$gp=paste(dso$GENEID, dso$dist_from_CDS_end, sep=':')

    #oligo.stats$gp=paste(oligo.stats$GENEID, oligo.stats$dist_from_CDS_end, sep=':')
    #oligo.stats$hmm2=dso$hmm2[match(oligo.stats$gp, dso$gp)]

    # integrate data from hmm 3 state into oligo table
    #so3=sapply(hmm.out, function(x) x$states.out3)
    #dso3=data.frame(GENEID=rep(names(so3), sapply(so3, length )), 
    #               dist_from_CDS_end=as.vector(unlist(sapply(so3, names))),
    #               hmm3=as.vector(unlist(so3)))
    #dso3$gp=paste(dso3$GENEID, dso3$dist_from_CDS_end, sep=':')
    #oligo.stats$hmm3=dso3$hmm3[match(oligo.stats$gp, dso3$gp)]
    #save(oligo.stats, file=paste0(out.base.dir, 'processed/RData/oligo.stats.2.RData'))
    #load(paste0(out.base.dir, 'processed/RData/oligo.stats.2.RData'))

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
formula.input=as.formula(slope.binarized~scale(cnt)+scale(p_intercept)+as.factor(expt)+scale(guide.GCcontent)+scale(PAMvariantCNT)+binot+U6.Terminator+Score+(1|GENEID)+(1|oligo))
oligo.stats=doGLMER.oligo.gene(formula.input, data.in=big.mm[!big.mm$drop,], oligo.stats, lab='EssentialAll')

formula.input=as.formula(slope.binarized~scale(cnt)+scale(p_intercept)+scale(guide.GCcontent)+scale(PAMvariantCNT)+binot+U6.Terminator+Score+(1|GENEID)+(1|oligo))
oligo.stats=doGLMER.oligo.gene(formula.input, data.in=big.mm[!big.mm$drop & (big.mm$expt == 'WT'),], oligo.stats, lab='EssentialWT')

formula.input=as.formula(slope.binarized~scale(cnt)+scale(p_intercept)+scale(guide.GCcontent)+scale(PAMvariantCNT)+binot+U6.Terminator+Score+(1|GENEID)+(1|oligo))
oligo.stats=doGLMER.oligo.gene(formula.input, data.in=big.mm[!big.mm$drop & (big.mm$expt == 'NMD'),], oligo.stats, lab='EssentialNMD')

formula.input=as.formula(slope.binarized~scale(cnt)+scale(p_intercept)+as.factor(expt)+scale(guide.GCcontent)+scale(PAMvariantCNT)+binot+U6.Terminator+Score+(1|oligo))
oligo.stats=doGLMER.oligo.gene(formula.input, data.in=big.mm[!big.mm$drop,], oligo.stats, lab='EssentialAll_noGene', doGene=FALSE)
#============================================================================================
save(oligo.stats, file=paste0(out.base.dir, 'processed/RData/oligo.stats.EssentialBlups.RData'))
load(paste0(out.base.dir, 'processed/RData/oligo.stats.EssentialBlups.RData'))


###### NEW ###################################################3
#remove gene effect
ibmm=glmer(slope.binarized~
                cnt+
                as.factor(expt)+
                p_intercept+
                dist_from_CDS_end+
                CDS_length+
                PAMvariantCNT+
                binot+
                (1|oligo),
                 family=binomial(link="logit"), data=big.mm[!big.mm$drop,], verbose=T)
rbmm=ranef(ibmm, cond=TRUE)
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




# Figure 1 pilot testing ------------------------------------------------------------------------------------
    git.dir='/media/jbloom/d1/coupled_CRISPR/Experiments/Preliminary/Pilot_092015/'
    load(paste0(git.dir, 'summaryStats.RData'))
    #load(paste0(git.dir, 'sampleTable.RData'))

    summary.stats.table=t(do.call('rbind', summary.stats))
    pilot.norm=t(t(summary.stats.table)/colSums(summary.stats.table))

    #1 = matches haplotype no mismatches or indels in 50bp window of pam
    #2 = indel in 50 bp window around pam
    #3 = everything else
    #4 = WT
    pilot.drop=c(7,9:16, 17:19 , 22:25) #9:16, 23:25)  #replot 22:25
    pilot.norm.drop=pilot.norm[,-pilot.drop][,-c(9:10)]
    pilot.norm.drop=pilot.norm.drop[,ncol(pilot.norm.drop):1]

    pdf(file='~/Desktop/Figure1B_Sri.pdf', width=22, height=8)
    par(oma=c(1,20,1,1))
    bp=barplot(pilot.norm.drop, col=c('black','red', 'grey', 'white'), # density=c(100,55,25,0),
               horiz=TRUE, xlim=c(0,1), yaxt='n', 
               beside=T, xlab='fraction of reads', space=c(0,3), width=2 )#, main='Figure1B')
    axis(2, bp[2,], colnames(pilot.norm.drop), las=2)
    legend('bottomright',fill=c('black','red', 'grey', 'white'), #density= #rev(c(100,55,25,0)) , 
                                               (c('matches expected edit', 
                                                    'indel near edit',
                                                    'mismatch near edit', 
                                                    'unedited')))
    dev.off()
    # go back and check this 030317
    # 4A-18A_S1:pJS4:chrV:33123 mismatches 
    # 4N-18N_S2:pJS4:chrV .. really so many indels??
#-------------------------------------------------------------------------------------------------------------



# Figure 2 slopes (Dubious vs essential)---------------------------------------------------------------------
    pdf(file='~/Desktop/CC/Figure2.pdf', width=15, height=15)
    par(oma=c(1,1,1,1))
    with(big.mm, {
        plot(jitter(as.numeric(dubious)), slope, col='#00000013', xaxt='n', xlab='', ylab='barcode slope', main='Figure 2B')
        axis(1,at=c(0,1), lab=c('essential', 'dubious'))  
        abline(h=-0.025, col='red')   
        text(c(.5,.5),c((-0.025-.02),(-0.025+.02)), c('dead', 'alive')) })
    dev.off()
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

# fit binarized slope model  Figure or table 4
b3=big.mm[!big.mm$drop & big.mm$dist_from_CDS_end<151, ] # & (big.mm$expt == 'WT') , ]#& big.mm$dist_from_CDS_end<300,]
b3$expt=as.factor(b3$expt)
b3$domain.downstream=as.factor(b3$domain.downstream)
b3$oligo=as.factor(b3$oligo)
# scale(H0_w)+
#scale(CDS_length)+
#binot+
# subset of data out to 125
# 2 step ... fit full model... then remove terms that explain nothin 
#mrna or some other data set

splic.genes=unlist(strsplit(GO.out[[1]][1,]$genes.in.set.sorted, ' '))
spgf=b3$GENEID %in% splic.genes

slopeAB=glmer(slope.binarized~
                cnt+
                as.factor(expt)+
                p_intercept+
                guide.GCcontent+
                dist_from_CDS_end+
                CDS_length+
                PAMvariantCNT+
                binot+
                lcd+
                evolvability+
                mean.aa.perfect.conserved+
                viable_annotation+
                end.conservation+
                domain.downstream+spgf+
                (1|GENEID)+(1|oligo), data=b3, 
            family=binomial(link="logit"),
            control=glmerControl(optimizer="Nelder_Mead", optCtrl = list(maxfun = 1e6)) , verbose=T)
            
    #summary(slope5)
    at=(anova(slopeAB, test='F'))
    length(predict(slopeA))
    at[order(at[,1], decreasing=T),]
ssc=summary(slopeAB)$coefficients

#

#gAp=plot(gammALL$gam)
#gApD=plot(gammALL.dd$gam)
#gApnD=plot(gammALL.ndd$gam)

#dev.off()
#plot(gAp[[1]]$x, gAp[[1]]$fit, type='l', xlim=c(200,0))
#points(gAp[[1]]$x, gAp[[1]]$fit+gAp[[1]]$se, type='l', lty=2)
#points(gAp[[1]]$x, gAp[[1]]$fit-gAp[[1]]$se, type='l', lty=2)

# raw data and lines from regression model


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
















