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

# for GO enrichment
library(topGO)
library(org.Sc.sgd.db)
library(fdrtool)
library(sisus)

library("BSgenome.Scerevisiae.UCSC.sacCer3")
sacCer3=BSgenome.Scerevisiae.UCSC.sacCer3
library("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene")
txdb=TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
library(VariantAnnotation)

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
giant.table=readRDS(file = paste0(out.base.dir, 'processed/RData/giant.table.RData'))
rownames(oligos)=names(seq.tables)

#library(stringr)
# create 'vcf outut

# Effect of downstream codon on oligo essentiality 
#ref.revcomp=as.vector(sapply(oligos$REFCODON, revcomp))
#alt.revcomp=as.vector(sapply(oligos$VARCODON, revcomp))
#oligo.refout=ifelse(oligos$codingStrand=='+', oligos$REFCODON, ref.revcomp)
#oligo.altout=ifelse(oligos$codingStrand=='+', oligos$VARCODON, alt.revcomp)
#vcf.out=data.frame(CHROM=oligos$chr, POS=oligos$start,  ID=rownames(oligos), REF=oligo.refout, ALT=oligo.altout, QUAL='', FILTER='', INFO='', stringsAsFactors=F)
#vcf.out[order(vcf.out$CHROM, vcf.out$POS),]->vcf.out
#write.table(vcf.out, file=paste0('/media/jbloom/d1/coupled_CRISPR/Reference/stops_table.vcf'), row.names=F, quote=F, sep='\t')
#java -Xmx4g -jar snpEff.jar -v R64-1-1.86 /media/jbloom/d1/coupled_CRISPR/Reference/stops_table.vcf > test.vcf                   

# try to extract downstream codon
#oligo.revcomp=as.vector(sapply(oligos$stopPAM_oligos, revcomp))
#oligo.inframe=ifelse(oligos$codingStrand=='+', oligos$stopPAM_oligos, oligo.revcomp)
#head(oligo.inframe,20)
#oif.table=sapply(oligo.inframe, s2c)
#downstream.codon=apply(oif.table[53:55,],2, c2s)
#downstream.base=(oif.table[56,])

#downstream.codon.factor=downstream.codon[match(rownames(red.effs$oligo), rownames(oligos)) ]
#downstream.base.factor=downstream.base[match(rownames(red.effs$oligo), rownames(oligos)) ]

#stripchart(red.effs$oligo[,1]~downstream.codon.factor, vertical=T, pch=21, method='jitter')

#dcf.c=translate(sapply(downstream.codon.factor,s2c))
#dcf.c[grep('\\*', dcf.c)]='STOP'
#summary(lm(red.effs$oligo[,1]~downstream.codon.factor))
#summary(lm(red.effs$oligo[,1]~dcf.c))
#summary(lm(scale(red.effs$oligo[,1])~downstream.base.factor-1))
#summary(lm(scale(red.effs$oligo[,1])~downstream.base.factor-1))
#stripchart(red.effs$oligo[,1]~downstream.base.factor, vertical=T, pch=21, method='jitter', col='#00000022')
#boxplot(red.effs$oligo[,1]~downstream.base.factor)
#, vertical=T, pch=21, method='jitter', col='#00000022')

# slope blups and model in position effect



#giant.table.bad=giant.table[giant.table$cigar!='101M',]
#giant.table.bad.more=cbind(giant.table.bad, oligos[match(giant.table.bad$oligo, rownames(oligos)),])
#giant.table.bad.more.filt=giant.table.bad.more[rowSums(giant.table.bad.more[,c(9:13,15:19)])>35,]
# For Meru to investigate indel oligos
#saveRDS(giant.table.bad.more.filt, file='~/Dropbox/Public/giant.table.w.bad.oligos.RDS')

#[1] 215449
#
#gtbs=(giant.table.bad[which(!grepl('S', giant.table.bad$cigar)),])
#"73H48M"

# was 35
cutoff_t0=20
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
expt.counts.red2=colSums(giant.table[(giant.table$WT_0> cutoff_t0 | giant.table$NMDm_0> cutoff_t0) & giant.table$cigar=='101M'  ,expt.columns])

barplot(rbind(expt.counts.red2, expt.counts.red-expt.counts.red2, expt.counts-expt.counts.red),
        ylab='total reads', xlab='sample', main='reads per sample', col=c('black', 'grey', 'white'))
legend('topright', c( paste0(round(sum(expt.counts)/1e6,1), ' x 1e6 total reads'),
                      paste0(round(sum(expt.counts.red)/1e6,1), ' x 1e6 total reads (reads per barcode filter)')  ,
                      paste0(round(sum(expt.counts.red2)/1e6,1), ' x 1e6 total reads (reads per barcode and perfect match repair)')),
            fill=c('black', 'grey', 'white'),
               cex=1.5)  


giant.table.RPKM=10^9*t(t(giant.table[,expt.columns])/(101*expt.counts))
colnames(giant.table.RPKM)=paste0(colnames(giant.table.RPKM), '_rpkm')
giant.table=cbind(giant.table, giant.table.RPKM)

giant.table$matched.barcode=giant.table$WT_0>cutoff_t0 & giant.table$NMDm_0>cutoff_t0
#215449
#R> sum(giant.table$cigar=='101M')
#[1] 348288

giant.table.WT  = giant.table[giant.table$WT_0>cutoff_t0,]
giant.table.WT  = giant.table.WT[,-which(grepl('^NMD', colnames(giant.table.WT)))]
giant.table.NMD = giant.table[giant.table$NMDm_0>cutoff_t0,]
giant.table.NMD  = giant.table.NMD[,-which(grepl('^WT', colnames(giant.table.NMD)))]

# extract oligos that perfectly sequence that was synthesized
gt.WT.p   = giant.table.WT[giant.table.WT$cigar=='101M',]
gt.NMD.p = giant.table.NMD[giant.table.NMD$cigar=='101M',]

# for example, find oligos with indel near cutting site 
#gt.WT.del=(giant.table.WT[grep('4\\dM\\dD|5\\dM\\dD', giant.table.WT$cigar),])
#gt.WT.del=cbind(gt.WT.del, oligos[match(gt.WT.del$oligo, rownames(oligos)),])

gt.WT.pe=cbind(gt.WT.p, oligos[match(gt.WT.p$oligo, rownames(oligos)),])
gt.WT.pe=cbind(gt.WT.pe, abundance[match(gt.WT.pe$GENEID, abundance$yORF),])

gWg=split(gt.WT.pe, gt.WT.pe$GENEID)
gWgL=lapply(gWg, function(x) split(x, x$dist_from_CDS_end))

gt.NMD.pe=cbind(gt.NMD.p, oligos[match(gt.NMD.p$oligo, rownames(oligos)),])
gt.NMD.pe=cbind(gt.NMD.pe, abundance[match(gt.NMD.pe$GENEID, abundance$yORF),])

gNg=split(gt.NMD.pe, gt.NMD.pe$GENEID)
gNgL=lapply(gNg, function(x) split(x, x$dist_from_CDS_end))

#this might not generalize if we change t0 count threshold
remove.ind=-which(!(names(gNgL) %in% names(gWgL)))
if(length(remove.ind)>0) { gNgL=gNgL[remove.ind] }

# Exploratory
#plot(gt.WT.pe$dist_from_CDS_end, gt.WT.pe$WT_48_rpkm/gt.WT.pe$WT_0_rpkm, col='#00000022')
#plot(gt.WT.pe$dist_from_CDS_end, gt.WT.pe$WT_48_rpkm/gt.WT.pe$WT_0_rpkm, col='#00000011', ylim=c(0,2))
#md=sapply(split( gt.WT.pe$WT_48_rpkm/gt.WT.pe$WT_0_rpkm, gt.WT.pe$dist_from_CDS_end), mean, na.rm=T)
#md.ceiling=md
#md.ceiling[md.ceiling>2]=1.5+rnorm(length(md[md>2]), .05, .05)
#plot(as.numeric(names(md)), as.vector(md.ceiling), ylim=c(0,2), xlim=c(0,1000))
#plot(as.numeric(names(md)), as.vector(md.ceiling), ylim=c(0,2), xlim=c(0,1000))
#md2=sapply(split( gt.NMD.pe$NMDm_48_rpkm/gt.NMD.pe$NMDm_0_rpkm, gt.NMD.pe$dist_from_CDS_end), mean, na.rm=T)
#md2.ceiling=md2
#md2.ceiling[md2.ceiling>2]=1.5+rnorm(length(md2[md2>2]),.05, .05,)
#points(as.numeric(names(md2)), as.vector(md2.ceiling), ylim=c(0,2), col='red')
#gt=cbind(giant.table, oligos[match(giant.table$oligo, rownames(oligos)),])
#gt=gt[gt$cigar=='101M',]
#gtm=gt[gt$matched.barcode,]
#gss=split((gtm$WT_48_rpkm/gtm$WT_0_rpkm)-(gtm$NMDm_48_rpkm/gtm$NMDm_0_rpkm), gtm$dist_from_CDS_end)
#stripchart(gss, vertical=T, pch=21)
#gtms=sapply(split((gtm$WT_48_rpkm/gtm$WT_0_rpkm)-(gtm$NMDm_48_rpkm/gtm$NMDm_0_rpkm), gtm$dist_from_CDS_end), median, na.rm=T)
#test=lapply(gWgL, function(x)  lapply(x, function(y)  y$WT_48_rpkm/y$WT_0_rpkm     )
#test2=lapply(gNgL, function(x)  lapply(x, function(y)  y$WT_48_rpkm/y$WT_0_rpkm     )

# the anova is the slow part, most fits are significant, can leave that out 
calc_pois_regression=function(counts, t.var) {
    tot.counts=colSums(counts)
    #bad.exp=c(5:6) #which(tot.counts<2e6)
    #counts[,bad.exp]=NA
    #nm=mcMap(function(i){
    #                glm(counts[i,]~offset(log(tot.counts)), family='poisson')
    #            },   1:nrow(counts),mc.cores=30)
    fm=mcMap(function(i){
                    glm(counts[i,]~t.var+offset(log(tot.counts)), family='poisson')
                },   1:nrow(counts),mc.cores=50)
 
    #am=simplify2array(mcMap(function(i){
    #           anova(fm[[i]],nm[[i]], test='Chisq')$Deviance[2] 
    #            }, 1:nrow(counts), mc.cores=30))
    #return(list(nm=nm, fm=fm,am=am))
    return(list(fm=fm))       
}

WT.fit=calc_pois_regression(data.matrix(gt.WT.pe[,9:13]), t.var)
WT.coef=sapply(WT.fit$fm, function(x) coef(x))
#library(FNN)

# look at breakdown of poisson scores per gene 
NMD.fit=calc_pois_regression(data.matrix(gt.NMD.pe[,9:13]), t.var)
NMD.coef=sapply(NMD.fit$fm, function(x) coef(x))
#library(FNN)

# Coefficient crashes to small value if counts at t24 =0
cb.thresh=-0.025
par(mfrow=c(2,1))
hist(WT.coef[2,], breaks=250, xlim=c(-.2, .05), main='WT', xlab='barcode slope')
abline(v=cb.thresh, col='red')
hist(NMD.coef[2,], breaks=250, xlim=c(-.2, .05), main='NMD', xlab='barcode slope')
abline(v=cb.thresh, col='red')


# Different way of getting depletion classfication ---------------
    rpkm.columns=grep('rpkm', names(gt.WT.pe))[1:5]
    gWw=gt.WT.pe[,rpkm.columns]
    #
    ## depletion 
    gWdep=gWw[,3]/apply(gWw[,c(1,2)], 1, max)
    gWdep.comp=log2(gWdep+.5)
    gWdep.binary2=1-((gWdep>.1)+0)
    #
    #rpkm.columns=grep('rpkm', names(gt.NMD.pe))[1:5]
    gNw=gt.NMD.pe[,rpkm.columns]
    #
    ## depletion 
    gNdep=gNw[,3]/apply(gNw[,c(1,2)], 1, max)
    gNdep.comp=log2(gNdep+.5)
    gNdep.binary2=1-((gNdep>.1)+0)
#-----------------------------------------------------------------

##~ .8 correlation with regression based classifier
gWdep.binary=(WT.coef[2,]>(cb.thresh))+0
gNdep.binary=(NMD.coef[2,]>(cb.thresh))+0

##rename columns containing count and rpkm data
colnames(gt.WT.pe)=gsub('WT', 't', colnames(gt.WT.pe))
colnames(gt.NMD.pe)=gsub('NMDm', 't', colnames(gt.NMD.pe))

## both experiments combined
atog=rbind(gt.WT.pe, gt.NMD.pe) # gt.WT.pe[,-grep('WT', colnames(gt.WT.pe))], gt.NMD.pe[,-grep('NMD', colnames(gt.NMD.pe))])
atog$intercept=c(WT.coef[1,], NMD.coef[1,])
atog$p_coeff=c(WT.coef[2,], NMD.coef[2,])
atog$dep=c(gWdep.binary, gNdep.binary)
atog$dep2=c(gWdep.binary2, gNdep.binary2)
atog$expt=c(rep('WT', length(gWdep.binary)), rep('NMD', length(gNdep.binary)))


## model on full data set 
#fModel0=glmer(dep~1+(1|GENEID),family=binomial(link="logit"), data=atog)
#fModel1=glmer(dep~1+(1|GENEID)+(1|oligo)+(1|expt),family=binomial(link="logit"), data=atog)
#fModel2=glmer(dep~1+(1|GENEID)+(1|oligo)+(1|expt)+(1|dubious),family=binomial(link="logit"), data=atog)

#+(1|expt)
fModel0=glmer(dep~1+intercept+(1|GENEID),family=binomial(link="logit"), data=atog)
fModel1=glmer(dep~1+intercept+(1|oligo),family=binomial(link="logit"), data=atog)
fModel3=glmer(dep~1+intercept+(1|GENEID)+(1|oligo),family=binomial(link="logit"), data=atog)

#rfull0=ranef(fModel0,condVar=TRUE)
#rfull1=ranef(fModel1,condVar=TRUE)
#rfull2=ranef(fModel2, condVar=TRUE)
rfull3=ranef(fModel3, condVar=TRUE)

par(mfrow=c(2,1))
hist(rfull3$oligo[[1]], breaks=100, xlab='oligo blups', main='oligo blups')
eo=exp(as.vector(rfull3$oligo[,1]))
lgit=eo/(1+eo) #
#lgit=as.vector(1/(1+exp(-rfull3$oligo[[1]]))) #/(1+exp(rfull3$oligo[[1]])))

plot(as.vector(rfull3$oligo[,1]), lgit , xlab='blup', ylab='exp(blup)/(1+exp(blup))', type='p', ylim=c(0,1))

par(mfrow=c(2,1))
hist(rfull3$GENEID[[1]], breaks=100, xlab='gene blups', main='gene blups')
lgit=as.vector(exp(rfull3$GENEID[[1]])/(1+exp(rfull3$GENEID[[1]])))
plot(as.vector(rfull3$GENEID[[1]]), lgit , xlab='blup', ylab='exp(blup)/(1+exp(blup))', type='p', ylim=c(0,1))

cor.test(predict(fModel3,type='response'), (atog$dep))

gdub=oligos$dubious[match(rownames(rfull3$GENEID), oligos$GENEID)]
wilcox.test(rfull3$GENEID[,1]~gdub)
#stripchart(rfull3$GENEID[,1]~gdub, vertical=T)

m.all=merge(rfull3$oligo, oligos, by='row.names')

# effect of distance vs oligo essentiality for the non-essential genes  (NULL)
plot(m.all$dist_from_CDS_end[m.all$dubious], m.all[m.all$dubious,2], xlim=c(0,150), ylab='oligo blup', xlab='dist from CDS end')
boxplot(m.all[m.all$dubious & m.all$dist_from_CDS_end<151 ,2]~m.all$dist_from_CDS_end[m.all$dubious & m.all$dist_from_CDS_end<151],xlim=c(0,150), add=T, xaxs=F)
points(smooth.spline(m.all$dist_from_CDS_end[m.all$dubious], m.all[m.all$dubious,2]), type='l', col='red')

ofs=overlapping_features[sapply(overlapping_features, length)>1]
ofsd=do.call('rbind', lapply(ofs, function(x) x[c(1,2)]))
ofsd=unique(rbind(ofsd, ofsd[,c(2,1)]))

gdub.names=rownames(rfull3$GENEID)[gdub]
#plot(frac.aa.perfect.conserved[as.vector(ofsd[,2][match(gdub.names, ofsd[,1])])], (rfull3$GENEID[gdub,1]))

# gene blup for non-essential vs conservation of overlapping essential gene 
plot(frac.aa.perfect.conserved[as.vector(ofsd[,2][match(gdub.names, ofsd[,1])])], (rfull3$GENEID[gdub,1]), xlab='fraction perfectly conserved for overlapping essential gene', ylab='gene blup')
cor.test(frac.aa.perfect.conserved[as.vector(ofsd[,2][match(gdub.names, ofsd[,1])])], (rfull3$GENEID[gdub,1]), method='spearman' )
plot(dn_ds$H0_w[match(as.vector(ofsd[,2][match(gdub.names, ofsd[,1])]), dn_ds$Gene)], (rfull3$GENEID[gdub,1]) )
cor.test(dn_ds$H0_w[match(as.vector(ofsd[,2][match(gdub.names, ofsd[,1])]), dn_ds$Gene)], (rfull3$GENEID[gdub,1]), method='spearman' )

#par(mfrow=c(1,2))
#stripchart(rfull1$GENEID[,1]~gdub, vertical=T, xlab='DUBIOUS ORF', ylab='persisitence score (+ = persisting)', jitter=.1, method='jitter')
stripchart(rfull3$GENEID[,1]~gdub, vertical=T, xlab='DUBIOUS ORF', ylab='gene persisitence score (+ = persisting)', jitter=.1, method='jitter', pch=21)

# remove dubious or flagsyn
# and add in additional column of hand - annotated  non-essentials form Meru
# load genes.to.drop
load(file='/media/jbloom/d1/coupled_CRISPR/Reference/additional_genes_to_drop_not_essential.RData') 

genes.to.drop=as.vector(unlist(genes.to.drop))
atog$genes.to.drop=atog$GENEID %in% genes.to.drop

atog.red=atog[atog$dubious==FALSE & atog$flagsyn==FALSE & atog$genes.to.drop==FALSE,]

# recalculate oligo and gene effects
rModel=glmer(dep~1+scale(intercept)+(1|GENEID)+(1|oligo),family=binomial(link="logit"), 
             data=atog.red,glmerControl(optimizer='Nelder_Mead',optCtrl=list(maxfun=5e4)))
#lModel=lmer(scale(p_coeff)~1+intercept+(1|GENEID)+(1|oligo), data=atog.redm)
lModel=glmer(dep~1+intercept+(1|oligo:GENEID), family=binomial, data=atog.red)
#clmodel=cor.test(predict(lModel,type='response'), (atog.redm$p_coeff))
#,type='response'
#cRmodel=cor.test(predict(rModel,type='response'), (atog.red$dep))

#null for variance components
set.seed=100
gdep.random=replicate(100, {sample(atog.red$dep) }, simplify=FALSE)
gred.rvar=mclapply(gdep.random, function(y) {
                    z=atog.red
                    z$dep=y
                    gout= glmer(dep~1+intercept+(1|GENEID)+(1|oligo)+(1|expt),family=binomial(link="logit"), data=z)
                    red.effs=(ranef(gout, condVar=TRUE))
                    return(list(red.effs=red.effs, var.out=sapply(VarCorr(gout), function(x) attr(x, 'stddev'))^2)) 
                }, mc.cores=10)
vc.null=do.call('rbind', lapply(gred.rvar, function(x) x$var.out))

# normalize VC output
sf=(rvcresults[1]+rvcresults[2])/cRmodel$estimate^2 
boxplot(vc.null/sf, ylim=c(0,.5), ylab='variance explained', names=c('OLIGO', 'GENE', 'WT v NMD'))
rvcresults=sapply(VarCorr(rModel), function(x) attr(x, 'stddev'))^2
points(1,rvcresults[1]/sf, col='red', cex=3, pch="*")
points(2,rvcresults[2]/sf, col='red', cex=3, pch="*")
points(3,0, col='red', cex=3, pch="*")

red.effs=(ranef(rModel, condVar=TRUE))
#cor.test(predict(rModel,type='response'), (atog.red$dep))

# visualize oligo effect as a function of distance from end of gene
oreff=oligos[match(rownames(red.effs$oligo), rownames(oligos)),]
#par(mfrow=c(2,1))
# dn/ds table -------------------------------------
odf=data.frame(oblup=red.effs$oligo[,1], dist_from_end=oreff$dist_from_CDS_end, 
               end_conservation=end.conservation[match(rownames(red.effs$oligo), rownames(oligos))],
               oligos[match(rownames(red.effs$oligo), rownames(oligos)),], stringsAsFactors=F
               )
sodf=odf[odf$dist_from_end<200,]
anova(lm(sodf$oblup~sodf$dist_from_end+sodf$end_conservation))
#plot(sodf$end_conservation, sodf$oblup     )
#par(mfrow=c(2,1))
#plot(sodf$dist_from_end, sodf$end_conservation, xlim=c(0,200))
plot(sodf$dist_from_end, sodf$oblup, col=ifelse(sodf$end_conservation>.9, 'red', 'black'), xlim=c(0,150), cex=.5,
     ylab='oligo persistence', xaxt='n', xlab='distance in amino acids from the stop')
boxplot(sodf$oblup[sodf$end_conservation>.9]~sodf$dist_from_end[sodf$end_conservation>.9], col='#ff000022', border='#ff000011', add=T, xaxs=F)
boxplot(sodf$oblup[sodf$end_conservation<.9]~sodf$dist_from_end[sodf$end_conservation<.9], col='#00000011',  border='#00000011', add=T, xaxs=F)
nsodf=na.omit(sodf)
points(smooth.spline(nsodf$dist_from_end[nsodf$end_conservation>.9], nsodf$oblup[nsodf$end_conservation>.9]), type='l', col='red', lwd=2)
points(smooth.spline(nsodf$dist_from_end[nsodf$end_conservation<.9], nsodf$oblup[nsodf$end_conservation<.9]), type='l', col='black', lwd=2)

#plot(jitter(oreff$dist_from_CDS_end,.5), red.effs$oligo[,1], xlim=c(0,100), col='#00000088')


argdf=atog.red[match(rownames(red.effs$GENEID), atog.red$GENEID),]

# correlation between abundance and random effect of gene?
plot( log2(argdf$YEPD.mean),red.effs$GENEID[,1], ylab='gene persistence', xlab='log2(GFP)')
cr1=cor.test(log2(argdf$YEPD.mean), red.effs$GENEID[,1], method='spearman')
legend('topright', paste0('rho=', round(cr1$estimate,2), '      p< ', format.pval(cr1$p.value)))

# half life
plot(y=red.effs$GENEID[,1], x=log2(argdf$Corrected.Half.Life), ylab='gene persistence', xlab='log2( protein half-life)')
cr1=cor.test(log2(argdf$Corrected.Half.Life), red.effs$GENEID[,1], method='spearman')
legend('topright', paste0('rho=', round(cr1$estimate,2), '      p< ', format.pval(cr1$p.value)))

#conservation
plot(y=red.effs$GENEID[,1],x=mean.aa.perfect.conserved[rownames(red.effs$GENEID)], ylab='gene persistence', xlab='conservation score', xlim=c(0.6,1))
cr1=cor.test(mean.aa.perfect.conserved[rownames(red.effs$GENEID)], red.effs$GENEID[,1], method='spearman')
legend('topleft', paste0('rho=', round(cr1$estimate,2), '      p< ', format.pval(cr1$p.value)))

plot(y=red.effs$GENEID[,1],x=dn_ds$H0_w[match(rownames(red.effs$GENEID), dn_ds$Gene)], ylab='gene persistence', xlab='dN/dS', xlim=c(0,.35))
cr1=cor.test(dn_ds$H0_w[match(rownames(red.effs$GENEID), dn_ds$Gene)], red.effs$GENEID[,1], method='spearman')
legend('topleft', paste0('rho=', round(cr1$estimate,2), '      p< ', format.pval(cr1$p.value)))

#length
plot(y=red.effs$GENEID[,1],x=as.numeric(argdf$CDS_length), ylab='gene persistence', xlab='CDS length', xlim=c(0,2000))
cr1=cor.test(as.numeric(argdf$CDS_length), red.effs$GENEID[,1], method='spearman')
legend('topright', paste0('rho=', round(cr1$estimate,2), '      p< ', format.pval(cr1$p.value)))

#haploinsufficiency
mhaplo=merge(haploinf, red.effs$GENEID, by.x='orf', by.y='row.names')
plot(as.factor(mhaplo$slow_ypd_het), mhaplo[,19])
wilcox.test(mhaplo[,19]~as.factor(mhaplo$slow_ypd_het))
cor.test(mhaplo[,19],mhaplo$HET_AV, method='spearman')

plot(mhaplo$HET_AV ,mhaplo[,19], ylab='gene persistence', xlab='haploinsufficiency') #, method='spearman')
cr1=cor.test(mhaplo[,19],mhaplo$HET_AV, method='spearman')
legend('topleft', paste0('rho=', round(cr1$estimate,2), '      p< ', format.pval(cr1$p.value)))


via_filtered=viable_table
#via_filtered=via_filtered[via_filtered$Gene.Systematic.Name %in% rownames(red.effs$GENEID),]
via_merged=merge(viable_table, red.effs$GENEID, by.x='Gene.Systematic.Name', by.y='row.names', all.y=T)
via_merged$Phenotype[is.na(via_merged$Phenotype)]='inviable'

stripchart(via_merged[,11]~via_merged$Phenotype, vertical=TRUE, pch=21,method='jitter', xlim=c(0.5,2.5), ylab='gene persistence', xlab='any SGD annotation for viable',
           sub=paste0('wilcox p< ', format.pval(wilcox.test(via_merged[,11]~via_merged$Phenotype)$p.value)))
boxplot(via_merged[,11]~via_merged$Phenotype,add=T,xaxt='n')

wilcox.test(via_merged[,11]~via_merged$Phenotype)
#via_filtered=viable_table[viable_table$Experiment.Type!='systematic mutation set',]
#via_filtered=viable_table[viable_table$Experiment.Type=='classical genetics',]
svf=split(via_filtered[,-2], via_filtered$Gene.Systematic.Name)
svf=lapply(svf, function(x) apply(x, 1, paste, collapse=':'))
svf=sapply(svf, paste, collapse=';;')
cviable=rownames(red.effs$GENEID) %in% via_filtered$Gene.Systematic.Name
cviable=rownames(red.effs$GENEID) %in% via_filtered$Gene.Systematic.Name & via_filtered$Strain.Background=="Sigma1278b"

wilcox.test(red.effs$GENEID[,1]~cviable) 
identify(cviable, red.effs$GENEID[,1], rownames(red.effs$GENEID))

# negative correlation between gene conservation and gene blup
cor.test(frac.aa.perfect.conserved[rownames(red.effs$GENEID)], red.effs$GENEID[,1], method='spearman')
cor.test(frac.aa.perfect.conserved[rownames(red.effs$GENEID)], red.effs$GENEID[,1], method='spearman')
cor.test(mean.aa.perfect.conserved[rownames(red.effs$GENEID)], red.effs$GENEID[,1], method='spearman')

# dn/ds vs gene blup (low dn/ds = conserved)
#cor.test(dn_ds$H0_w[match(rownames(red.effs$GENEID), dn_ds$Gene)],  red.effs$GENEID[,1], method='spearman')
#cor.test(dn_ds$dN_H0[match(rownames(red.effs$GENEID), dn_ds$Gene)]/dn_ds$Length[match(rownames(red.effs$GENEID), dn_ds$Gene)],  
#         red.effs$GENEID[,1], method='spearman')
#cor.test(dn_ds$dN_H0[match(rownames(red.effs$GENEID), dn_ds$Gene)],  
#         red.effs$GENEID[,1], method='spearman')
#cor.test(dn_ds$dS_H0[match(rownames(red.effs$GENEID), dn_ds$Gene)],  
#         red.effs$GENEID[,1], method='spearman')
#cor.test(dn_ds$H0_w[match(rownames(red.effs$GENEID), dn_ds$Gene)],  red.effs$GENEID[,1], method='spearman')


# evolvability 
evolv.vec=evolvability$X[match(rownames(red.effs$GENEID), evolvability$Systematic.name)]
lm(red.effs$GENEID[,1]~evolv.vec)
stripchart(red.effs$GENEID[,1]~evolv.vec, vertical=T, method='jitter', pch=21, xlim=c(0.5,4.5), ylab='gene persistence')
boxplot(red.effs$GENEID[,1]~evolv.vec, add=T, xaxt='n') 
pairwise.t.test(red.effs$GENEID[,1],evolv.vec, p.adjust='bonferroni')      

# human complements null
ch.vec=comp_human$Final.CompStatus[match(rownames(red.effs$GENEID), comp_human$ScENSP)]
lm(red.effs$GENEID[,1]~ch.vec-1)
stripchart(red.effs$GENEID[,1]~ch.vec-1, vertical=T, method='jitter', pch=21, xlim=c(0.5,2.5), ylab='gene persistence', xlab='human CDS complements null')
boxplot(red.effs$GENEID[,1]~ch.vec, add=T, xaxt='n')

# difference based on AA changed to stop
mor=merge(oligos, red.effs$oligo, by='row.names')
mor1=mor[mor$dist_from_CDS_end<50,]
mors=split(mor, mor$guide)

stripchart(mor1$'(Intercept)'~mor1$REFAA,vertical=T, method='jitter',pch=21,ylab='oligo persistence', xlab='reference AA (edited to STOP)', main='edits within 50 AA of CDS end')
boxplot(mor1$'(Intercept)'~mor1$REFAA,add=T, xaxt='n', skinny=T, border='#00000077')

mor=merge(oligos, rfull3$oligo, by='row.names')
mors=split(mor, mor$guide)

mors2=mors[which(sapply(mors, function(x) nrow(x)>1))]
#mors2=mors2[-grep('AACACCGCTTCAAATATGCT', names(mors2))]
#mors2=mors2[-grep('ATTACCCATTTGAGATGCAT', names(mors2))]

dub.pers=as.vector(sapply(mors2, function(x) x$'(Intercept)'[x$dubious][1]))
ndub.pers=as.vector(sapply(mors2, function(x) x$'(Intercept)'[!x$dubious][1]))

group=as.factor(c(rep('dub', 53), rep('ndub', 53)))

plot((group), as.numeric(c(dub.pers, ndub.pers)))

points(as.numeric(group), as.numeric(c(dub.pers, ndub.pers)))
segments(1, dub.pers, 2, ndub.pers)
as.numeric(c(dub.pers, ndub.pers))
obg=split(oligos, oligos$guide)
osub=(oligos[oligos$guide %in% names(which(sapply(obg, function(x) nrow(x))>1)),])   
test=merge(osub, red.effs$oligo, by='row.names')


testGenes=red.effs$GENEID[,1] 
names(testGenes)=rownames(red.effs$GENEID)
        
geneSel.fx=function(allScore) { return(allScore<1e9) }
GO.out=list()
for(thisOntology in c('BP', 'MF', 'CC') ) {
        print(thisOntology)
        #thisOntology='BP'
        #thisOntology='MF'
        #thisOntology='CC'
        #GOData = new("topGOdata", ontology=thisOntology, allGenes = testGenes, annot = annFUN.gene2GO, gene2GO = gene2GOList, nodeSize=5)
        GOData = new("topGOdata", ontology=thisOntology, allGenes = testGenes, geneSel=geneSel.fx, annot = annFUN.gene2GO, gene2GO = gene2GOList, nodeSize=5)
        GOresult.dec = runTest(GOData, algorithm="classic", statistic="ks", scoreOrder='decreasing')
        GOresult.inc = runTest(GOData, algorithm="classic", statistic="ks", scoreOrder='increasing')

        gt=GenTable(GOData, GOresult.dec, numChar=140, topNodes = length(score(GOresult.dec)))
        gt2=GenTable(GOData, GOresult.inc, numChar=140, topNodes = length(score(GOresult.inc)))
        gt3=merge(gt,gt2,by='GO.ID')
        min.p=ifelse(as.numeric(gt3$result1.x)<as.numeric(gt3$result1.y), as.numeric(gt3$result1.x), as.numeric(gt3$result1.y))
        effect.direction=ifelse(as.numeric(gt3$result1.x)<as.numeric(gt3$result1.y), 'persisting', 'dead')
        
        gt3$min.p=as.numeric(min.p)
        gt3$effect.direction=effect.direction

        gt3=gt3[order(gt3$min.p, decreasing=F),]
        gt3=gt3[,c(1,2,3,12,13)]
        names(gt3)=c('GO.ID', 'Term', 'Annotated', 'p.value', 'class')

        GO.out[[thisOntology]]=gt3
        print(head(gt,20))
        print(head(gt2,20))
 }
library(qvalue)

fdr.out=qvalue(unlist(sapply(GO.out, function(x) x$p.value)), lambda=seq(.3,.5,.05), fdr.level=.01)
max(fdr.out$pvalues[fdr.out$significant])
#0.00086
thisOntology='CC'
GOData = new("topGOdata", ontology=thisOntology, allGenes = testGenes, geneSel=geneSel.fx, annot = annFUN.gene2GO, gene2GO = gene2GOList, nodeSize=5)
goID='GO:0005681'
plot(showGroupDensity(GOData, goID, ranks=T))
title('spliceosomal complex' )

thisOntology='MF'
GOData = new("topGOdata", ontology=thisOntology, allGenes = testGenes, geneSel=geneSel.fx, annot = annFUN.gene2GO, gene2GO = gene2GOList, nodeSize=5)
goID='GO:0003824'
plot(showGroupDensity(GOData, goID, ranks=T))
title('catalytic activity' )


 #       for ( i in 1:50) {
 #       goID=allRes[i,'GO.ID']
 #       goID='GO:0000375'
 #       plot(showGroupDensity(GOData, goID, ranks=T) )
 #       readline()
 #       }
 #       genesinterms=genesInTerm(GOData, goID)
 #       genes.enriched.list=lapply(genesinterms, function(x) x[x%in%names(testGenes[testGenes==1])])
 #       genes.enriched.list.simple=lapply(genes.enriched.list, function(x) as.character(SYS2ORF.key[x]))
 #       gt$Genes=as.vector(sapply(genes.enriched.list.simple, paste, collapse=','))
 #       gt$GenesSystematic=   as.vector(sapply(genes.enriched.list, paste, collapse=','))



atog.redm=atog.red
atog.redm$p_coeff[atog.red$p_coeff<(-.2)]=-.2
#library(depmixS4)
#gred=atog.redm[atog.redm$GENEID=='YER112W',]
#gred=atog.redm[atog.redm$GENEID=='YMR240C',]

#library(hmm.discnp)
#library(seqHMM)
#plot(gred$dist_from_CDS_end, gred$p_coeff)

#obs=split(gred$dep, gred$dist_from_CDS_end)

# set up matrices 
# get some prior probabilities for the emission matrix
#tot.E.oligos=length(atog$dep[atog$dubious | atog$flagsyn])
#(tot.E.oligos-sum(atog$dep[atog$dubious | atog$flagsyn]))/tot.E.oligos

#(alive)/total = (probability if alive seeing alive = 60.1%)
(sum(atog$dep[atog$dubious | atog$flagsyn]))/length(atog$dep[atog$dubious | atog$flagsyn])

tot.E.oligos2=length(atog$dep[!atog$dubious | !atog$flagsyn & (atog$dist_from_CDS_end>100)])

#  (probability if dead, seeing dead)  = 74%
(tot.E.oligos2-sum(atog$dep[!atog$dubious | !atog$flagsyn & (atog$dist_from_CDS_end>100)]))/tot.E.oligos2

# alive and really alive 
g.alive=c("YDR160W", "YKL141W", "YHR128W" ,"YJL156C",  "YKL037W", "YIL020C", "YFR029W", "YJR012C")
table(atog[atog$GENEID %in% g.alive, 'dep'])


# get prior probs for oligo closest to end --------------------------------------------------------------------------------------
#ogs=split(oligos, oligos$GENEID)
#og.closest=sapply(ogs, function(x) x$unique.Index[which.min(x$dist_from_CDS_end)])
#oblup.ind=rownames(red.effs$oligo)
#match(oblup.ind, rownames(oligos))
#R> sum(red.effs$oligo[na.omit(match(og.closest, match(oblup.ind, rownames(oligos)))),]>0)
#[1] 570
#R> length(red.effs$oligo[na.omit(match(og.closest, match(oblup.ind, rownames(oligos)))),]>0)
#[1] 909
##570/909
#---------------------------------------------------------------------------------------------------------------------------
#initialize HMM
# was .8 and .2
initprobs2=c(.63,.37)
transitionProbs2=rbind(c(.1,.9),
                       c(0,  1))
# top row was 60% and 40%
# from g.alive reset to 74% and 26%
emissionProbs2=rbind(c(.74, .26), 
                     c(.26, .74))

initprobs3=c(.63, .37, 0)
                         #A     D     #E
transitionProbs3=rbind( c(.1,   .9,      0),
                        c( 0,   .98,   .02),
                        c(0,     1,     0))

# rows are states, columns are observations
# two observable states and 3 hidden 
emissionProbs3=rbind(   c(.74,   .26, 0),
                       c(.26,  .74,  0 ),
                       c(.74,   .26, 0 ))
                    
atog.redm=atog.red
atog.redm$p_coeff[atog.red$p_coeff<(-.2)]=-.2
atog.redm$DA=as.factor(ifelse(atog.red$dep, '1', '2'))

atog.by.gene=split(atog.redm, atog.redm$GENEID)
atog.by.gene=lapply(atog.by.gene, function(x) {
                        x=x[order(x$dist_from_CDS_end),];
                        x$tt=as.numeric(as.factor(x$dist_from_CDS_end)); 
                        return(x);
            })

#gene='YMR240C'
#gene='YAL003W'

# print(gene)
#    abg=atog.by.gene[[gene]]#
#dm=depmix(dep~1, data=abg, nstates=2, family=binomial(), homogeneous=FALSE , trDens=transitionProbs, transition=~tt )
#fdm=fit(dm)
#depmixS4::posterior(fdm)
# mout=msm(DA~tt, data=abg,
#                     qmatrix=transitionProbs, ematrix=emissionProbs, exacttimes=T,
#                     obstype=rep(1, nrow(abg)), obstrue=rep(FALSE, nrow(abg)), initprobs=c(.5,.5) )
#viterbi.msm(mout)

#sdt=split(abg$DA, abg$tt)
#ymat=matrix(NA,nrow(abg),10)
#for(i in 1:nrow(abg)){
#    ymat[i, abg$tt[i]]=as.character(abg$DA[i])
#}
#hdo=hmm.discnp::hmm(ymat, yval=c('1','2'), K=2) # par0=list(tpm=transitionProbs, Rho=emissionProbs))
#hmm.discnp::viterbi(hdo)


doHMM.2state=function(abg, transitionProbs, emissionProbs, initprobs) {
    return(
        try(msm(DA~tt, data=abg,
                     qmatrix=transitionProbs,
                     hmodel=list(hmmCat(prob=c(emissionProbs[1,])), 
                                 hmmCat(prob=c(emissionProbs[2,]))), 
                     exacttimes=T, obstype=rep(1, nrow(abg)), 
                     obstrue=rep(FALSE, nrow(abg)), initprobs=initprobs) )
    )

}

doHMM.3state=function(abg, transitionProbs, emissionProbs, initprobs) {
    return(
        try(msm(DA~tt, data=abg,
                     qmatrix=transitionProbs,
                     hmodel=list(hmmCat(prob=c(emissionProbs[1,])), 
                                 hmmCat(prob=c(emissionProbs[2,])), 
                                 hmmCat(prob=c(emissionProbs[3,])) ),
                     exacttimes=T, obstype=rep(1, nrow(abg)), 
                     obstrue=rep(FALSE, nrow(abg)), initprobs=initprobs) )
    )
}


doBIC=function(l, k, n) { -2*l+k*log(n) }
#V5
#hmm.out=list()
hmm.out=foreach(gene = names(atog.by.gene) ) %do% {
    print(gene)
    abg=atog.by.gene[[gene]]# 'YMR240C']]
  
    # for null reverse the order of gene 
    abg.r=abg[order(abg$tt, decreasing=T),]
    abg.rle=rle(abg.r$dist_from_CDS_end)
    abg.r$tt=rep(1:length(abg.rle$lengths), abg.rle$lengths)
  
    initprobs2.new=initprobs2
    mout2=doHMM.2state(abg,transitionProbs2,emissionProbs2, initprobs2)
    while(class(mout2)!='msm') {
        ni=runif(1)
        initprobs2.new=c(ni,1-ni)
        mout2=doHMM.2state(abg,transitionProbs2,emissionProbs2, initprobs2.new)
    }
    mout2.rev = doHMM.2state(abg.r,transitionProbs2,emissionProbs2, initprobs2.new)
    
    mout3     = doHMM.3state(abg,transitionProbs3, emissionProbs3, initprobs3)
    #mout3.rev = doHMM.3state(abg.r,transitionProbs3, emissionProbs3, initprobs3)

    vout=viterbi.msm(mout2)
    vout3=NULL
    states.out3=NULL
    if(class(mout3)=='msm') {
        vout3=viterbi.msm(mout3)
        #dpoint3=abg$dist_from_CDS_end[min(which(vout$fitted==2))]
        x3=abg$dist_from_CDS_end
        states.out3=sapply(split(vout3$fitted, abg$dist_from_CDS_end), function(i) i[1])
    }
    dpoint=abg$dist_from_CDS_end[min(which(vout$fitted==2))]
    x=(cbind(vout, abg$dist_from_CDS_end))
    states.out=sapply(split(x$fitted, x[,6]), function(x) x[1])

    #hmm.out[[gene]]=
        return(list(mout=mout2, mnull=mout2.rev, vout=vout, mout3=mout3, vout3=vout3, states.out=states.out, states.out3=states.out3))
}
names(hmm.out)=names(atog.by.gene)
#save(hmm.out, file=paste0(out.base.dir, 'processed/RData/hmm.out.RData'))


# contrastbic
bic32=lapply(hmm.out, function(x) {
    if(class(x$mout)=='msm' & class(x$mout3)=='msm') {
        return(  list(bic3=as.numeric(doBIC(logLik.msm(x$mout3), nrow(x$mout3$data$mf), 6)),
                      bic2=as.numeric(doBIC(logLik.msm(x$mout), nrow(x$mout$data$mf), 3))) )
   } else{return(list(bic3=NULL, bic2=NULL))} })
b3=unlist(sapply(bic32, function(x)x$bic3))
 
b2=unlist(sapply(bic32, function(x)x$bic2))

# viterbi vs data
#

r.hmm=sapply(hmm.out, function(x) cor(x$vout$observed, x$vout$fitted))

hmm.scores=sapply(hmm.out, function(x) {
                      if(class(x$mout)=='msm' & class(x$mnull)=='msm') {
                        return(  logLik.msm(x$mout)-  logLik.msm(x$mnull) )
                      } else {return(NA) }
    } )



hmm.state.lods=sapply(hmm.out, function(x) {
    p.state=x$vout[,5]
    p.lod=log10(p.state[,1]+.0001)-log10(p.state[,2]+.0001)
    spl=split(p.lod, x$vout[,2])
    sapply(spl, function(x) x[length(x)])
    })

hmm.state.lods3=sapply(hmm.out, function(x) {
    p.state=x$vout3[,5]
    p.lod1=log10(p.state[,3]+.0001)-log10(p.state[,2]+.0001)
    p.lod2=log10(p.state[,3]+.0001)-log10(p.state[,1]+.0001)
    p.lod=ifelse(p.lod1>p.lod2, p.lod1, p.lod2)
    spl=split(p.lod, x$vout[,2])
    sapply(spl, function(x) x[length(x)])
    })


#df.hmm=data.frame(hmm.scores=hmm.scores, hmm.r=r.hmm, fraction.oligos.dead=sos)
#write.table(df.hmm, file='/media/jbloom/d1/coupled_CRISPR/output/position_dependent_effects_HMM4.txt', sep='\t' ,quote=F)

# important for output 
so=sapply(hmm.out, function(x) x$states.out)
so3=sapply(hmm.out, function(x) x$states.out3)

sapply(so, function(x) names(x)[x=='1'])
alive.pos=sapply(so, function(x) as.numeric(names(x)[x=='1']))
dead.pos=sapply(so, function(x) as.numeric(names(x)[x=='2']))
sos=sapply(so, function(x) sum(x-1)/length(x))


#976 where it could compute anything
hist(hmm.scores[sos>0 & sos<.99])

write.table(sort(hmm.scores[sos>0 & sos<.99], decreasing=T), file='~/a.txt', sep='\t' ,quote=F)

alive.cnt=sapply(alive.pos,length)
dead.cnt=sapply(dead.pos,length)

mean.aa.perfect.conserved[names(mean.aa.perfect.conserved) %in% names(alive.cnt)]
cor.test(alive.cnt, mean.aa.perfect.conserved[names(mean.aa.perfect.conserved) %in% names(alive.cnt)] ,method='spearman')

mean.aa.perfect.conserved[names(mean.aa.perfect.conserved) %in% names(dead.cnt)]
cor.test(alive.cnt/dead.cnt, mean.aa.perfect.conserved[names(mean.aa.perfect.conserved) %in% names(dead.cnt)] ,method='spearman')

plot(alive.cnt/dead.cnt, mean.aa.perfect.conserved[names(mean.aa.perfect.conserved) %in% names(dead.cnt)])

red.effs$oligo

i.overlap=list()
for(g in names(so) ){
    g.l=oligos[match(g, oligos$GENEID),'CDS_length']
    # positions called alive
    da1=as.numeric(names(so[[g]])) #which(so[[g]]==1)))
    da=g.l-da1
    o.blups.ind=grep(paste0(':',g,':'), rownames(red.effs$oligo))
    o.blups=red.effs$oligo[o.blups.ind,]
    o.blup.names=rownames(red.effs$oligo)[o.blups.ind]
    o.blup.order=as.numeric(sapply(strsplit(o.blup.names, ':'), function(x) x[24] ))

    o.blups=o.blups[match(da1, o.blup.order)]
    #if(sum(onames)!=length(da) ) {
    #    print(paste(g, 'error'))
    #}
    
    hmm.lods=as.vector(hmm.state.lods[[g]])
    hmm.lods3=as.numeric(as.vector(hmm.state.lods3[[g]]))
    so3.out=so3[[g]]
    if(is.null(so3.out)) {so3.out=hmm.lods3 }
    #domains
    if(!is.null(dgsplit[[g]])) {
        dg=Intervals(dgsplit[[g]][,c(6,7)])
        io=interval_overlap(da, dg)
        d.overlap=sapply(io, function(x) dgsplit[[g]][x,4])
        d.overlap.desc=sapply(io, function(x) dgsplit[[g]][x,5])
        d.overlap.start =sapply(io, function(x) dgsplit[[g]][x,6])
        d.overlap.stop = sapply(io, function(x) dgsplit[[g]][x,7])

        d.overlap=sapply(d.overlap, function(x) paste0(x, collapse=':'))
        d.overlap.desc=sapply(d.overlap.desc, function(x) paste0(x, collapse=':'))
        d.overlap.start=sapply(d.overlap.start, function(x) paste0(x, collapse=':'))
        d.overlap.stop=sapply(d.overlap.stop, function(x) paste0(x, collapse=':'))

        d.table=data.frame(gene=rep(g,length(da)), dist_from_end=da1, dist_from_start=da, oligo.blup=o.blups, hmm.call=so[[g]], hmm.lods=hmm.lods,
                           hmm.call3=so3.out, hmm.lods3=hmm.lods3,
                           domain=d.overlap, domain.desc=d.overlap.desc, 
                           domain.start=d.overlap.start, domain.stop=d.overlap.stop,
                           stringsAsFactors=F)
        i.overlap[[g]]=d.table
    }
    else{
        d.table=data.frame(gene=rep(g,length(da)), dist_from_end=da1, dist_from_start=da, oligo.blup=o.blups,  hmm.call=so[[g]], hmm.lods=hmm.lods,
                           hmm.call3=so3.out, hmm.lods3=hmm.lods3,
                           domain='',domain.desc='',
                           domain.start=NA, domain.stop=NA,
                           stringsAsFactors=F)
        i.overlap[[g]]=d.table
    }

}

i.df=do.call('rbind', i.overlap)
i.df.s=split(i.df, i.df$domain)
sort(sapply(i.df.s, function(x) sum(x$hmm.call==1)/length(x$hmm.call)))
i.df.g=split(i.df, i.df$gene)

i.df.g=lapply(i.df.g, function(x) {
       ii=x
       ii$domain.downstream=cumsum(ii$domain!='')>0
       return(ii)
})
i.df.g.df=do.call('rbind', i.df.g)
write.table(i.df.g.df, file='/media/jbloom/d1/coupled_CRISPR/output/oligo_blups_and_hmm_v5.txt',  quote=FALSE, sep='\t', row.names=FALSE)

t.test(i.df.g.df$hmm.call[i.df.g.df$domain.downstream], i.df.g.df$hmm.call[!i.df.g.df$domain.downstream])
boxplot(i.df.g.df$hmm.call~i.df.g.df$domain.downstream)
stripchart(i.df.g.df$hmm.call~i.df.g.df$domain.downstream, vertical=T, method='jitter')
plot(jitter(i.df.g.df$hmm.call[i.df.g.df$domain.downstream]),i.df.g.df$hmm.lods[i.df.g.df$domain.downstream])

#sum(ubl>3 & ubl <(-3))
#[1] 0
#R> sum(ubl<3 & ubl >(-3))
#[1] 1952
#R> length(ubl)
#[1] 8257
#R> 1952/8257
#[1] 0.2364


no.domain.genes=names(which(sapply(i.df.g, function(x) sum(x$domain!=""))==0))

wilcox.test(red.effs$GENEID[rownames(red.effs$GENEID) %in% no.domain.genes,], red.effs$GENEID[!(rownames(red.effs$GENEID) %in% no.domain.genes),])



t.var=c(0,24,48,72,96)
#pdf(file='~/Desktop/plots.pdf', width=10,height=12)
dir.create('/media/jbloom/d1/coupled_CRISPR/plots/HMM_v5/')
for(i in 1:length(gWgL)) {
    gene=names(gWgL)[i]
    dub.stat=ifelse(oligos[oligos$GENEID==names(gWgL)[i],]$dubious[1], 'DUBIOUS', 'ESSENTIAL')
    pdf(file=paste0('/media/jbloom/d1/coupled_CRISPR/plots/HMM_v5/', gene, ':', dub.stat, '.pdf'), width=12, height=13)

    dub.stat=ifelse(oligos[oligos$GENEID==names(gWgL)[i],]$dubious[1], 'DUBIOUS', 'ESSENTIAL')
    gene=names(gWgL)[i]
    print(i)
  
    conserv.vector=conservation[[gene]]

    cds.length=gWgL[[i]][[1]]$CDS_length[1]
    wt.length=length(gWgL[[i]])
    nmd.length=length(gNgL[[i]])

    wt.pos=as.numeric(names(gWgL[[i]]))
    nmd.pos=as.numeric(names(gNgL[[i]]))

    noligos=unique(sort(c(nmd.pos, wt.pos)))
    total.oligos=length(noligos)

    m=cbind(c(2:(total.oligos+1)),c((total.oligos+2):(total.oligos+total.oligos+1)))
    m=rbind(c(1,1), m)
    layout(m)
    par(mar=c(3,1,1,1), oma=c(2,1,3,1))
    if( is.null(conserv.vector) ) {  # length(coding.file)==0 | is.na(coding.file)){
        plot(0,0, xlim=c(1,cds.length), ylim=c(0,1), type='n') 
    } else {  
           # r=bio3d::read.fasta(coding.file)
           # tryCatch( {  plot.fasta(r, scale=3, aln.col='lightgrey', row.spacing=.2, plot.axis=F)}, error=function(e) {} )
            #plot(0,0, xlim=c(1,cds.length), ylim=c(0,8), type='n') 
           plot(conserv.vector, type='l', xlim=c(1,cds.length), ylim=c(0,1), xaxt='n', ylab='fraction aa conserved', col='#00000066', lwd=2)
           #lim <- par("usr")
           #rect(lim[1]-1, lim[3]-1, lim[2]+1, lim[4]+1, border = "red", col = "pink")
           #points(conserv.vector, type='h', xlim=c(1,cds.length), ylim=c(0,1), xaxt='n', ylab='fraction aa conserved')

    }
# color ticks by HMM output
    hmm.calls=so[[gene]]
    hmm.calls3=so3[[gene]]

    unique.xaxs.locations=unique(c(wt.pos, nmd.pos))
   # text(cds.length-wt.pos, .5, wt.pos, srt=90, col='red', cex=1.5)
   # text(cds.length-nmd.pos, .5, nmd.pos, srt=90, col='red', cex=1.5)
    axis(1, at=c(cds.length-unique.xaxs.locations), labels = unique.xaxs.locations, cex.axis=1) 
     #at=c(1, cds.length))
    if(!is.null(dgsplit[[gene]])) {
        rect(dgsplit[[gene]][,6], 0, dgsplit[[gene]][,7], 1, col='#0000ff33')
        text(dgsplit[[gene]][,6]+20,.3, dgsplit[[gene]][,5], cex=1.2)
    }
   if(!is.null(hmm.calls)){
        if(length(which(hmm.calls==1))>0 ) { }
            yp=cds.length-as.numeric(names(hmm.calls)[which(hmm.calls==1)])
            points(yp, rep(.1, length(yp)), type='l', col='green', pch=20, cex=2, lwd=4, lty=1)
            points(yp, rep(.1, length(yp)), type='p', col='green', pch=20, cex=2, lwd=4, lty=1)

        if(length(which(hmm.calls==2))>0){
            yp=cds.length-as.numeric(names(hmm.calls)[which(hmm.calls==2)])
            points(yp, rep(.1, length(yp)), type='l', col='red', pch=20, cex=2, lwd=4, lty=1)
            points(yp, rep(.1, length(yp)), type='p', col='red', pch=20, cex=2, lwd=4, lty=1)

        }
    }


# extract conservation track get info for scer coords 
    #par(mfrow=c(length(gWgL[[i]]), 1), 
    #par(mar=c(2.5,4,1,0), oma=c(1,1,3,1))
    for(jj in 1:total.oligos) {
        if( noligos[[jj]] %in% wt.pos) {
            j=match(noligos[jj], wt.pos) 
            plot(t.var, gWgL[[i]][[j]][1,9:13], type='n', xlab='time', ylab='RPKM', ylim=c(0, max( gWgL[[i]][[j]][,9:13])), main=names(gWgL[[i]])[j] )
            lim <- par("usr")
            if(!is.na(match(noligos[jj], names(hmm.calls))) ){
               mind=match(noligos[jj], names(hmm.calls) )
                 if(hmm.calls[[mind]]==2) {
                     rect(lim[1]-1, lim[3]-1, lim[2]+1, lim[4]+1,  col = "#ff000011")
                  } 
                  if(hmm.calls[[mind]]==1) {
                      rect(lim[1]-1, lim[3]-1, lim[2]+1, lim[4]+1, col = "#00ff0011")
                 }
            }
             if(!is.na(match(noligos[jj], names(hmm.calls3))) ){
                 mind=match(noligos[jj], names(hmm.calls3) )
                 if(hmm.calls3[[mind]]==3) {
                 legend('topright', 'E', text.col='brown', cex=1.5, bty='n')
                }
             }
            if( sum(gWgL[[i]][[j]]$flagsyn)>0) { legend('topright', 'flagsyn', text.col='red') }
            for(k in 1:nrow(gWgL[[i]][[j]]) ){ points(t.var, gWgL[[i]][[j]][k,9:13], type='b', xlab='time', ylab='RPKM', col=gWgL[[i]][[j]][k,'matched.barcode']+1 )          }
        } else {plot(0,0, type='n', axes=FALSE) }
     }
    for(jj in 1:total.oligos) {
        if( noligos[[jj]] %in% nmd.pos) {
            j=match(noligos[jj], nmd.pos) 
            plot(t.var, gNgL[[i]][[j]][1,9:13], type='n', xlab='time', ylab='RPKM', ylim=c(0, max( gNgL[[i]][[j]][,9:13])), main=names(gNgL[[i]])[j] )
            lim <- par("usr")
             if(!is.na(match(noligos[jj], names(hmm.calls))) ){
               mind=match(noligos[jj], names(hmm.calls) )
                 if(hmm.calls[[mind]]==2) {
                     rect(lim[1]-1, lim[3]-1, lim[2]+1, lim[4]+1,  col = "#ff000011")
                  } 
                  if(hmm.calls[[mind]]==1) {
                      rect(lim[1]-1, lim[3]-1, lim[2]+1, lim[4]+1, col = "#00ff0011")
                 }
            }
             if(!is.na(match(noligos[jj], names(hmm.calls3))) ){
                 mind=match(noligos[jj], names(hmm.calls3) )
                 if(hmm.calls3[[mind]]==3) {
                 legend('topright', 'E', text.col='brown',cex=1.5, bty='n')
                     }
             }

            if( sum(gNgL[[i]][[j]]$flagsyn)>0) { legend('topright', 'flagsyn', text.col='red') }
            for(k in 1:nrow(gNgL[[i]][[j]]) ){ points(t.var, gNgL[[i]][[j]][k,9:13], type='b', xlab='time', ylab='RPKM', col=gNgL[[i]][[j]][k,'matched.barcode']+1)          }
        } else {plot(0,0, type='n', axes=FALSE) }

    }
    if(!is.na(match(gene, names(gene2alias))) ) { #null(gene2alias[[gene]])) {
        title(paste(gene, ',', dub.stat, '\n', gene2alias[[gene]]), outer=T)
    } else {
        title(paste(gene, ',', dub.stat), outer=T) 
    }


    dev.off()

}
























# load Provean Scores
load('/data/Databases/provean_scores.RData')

#load grantham matrix
grantham=read.delim('/media/jbloom/d1/coupled_CRISPR/Reference/grantham.matrix.csv', header=T, sep='\t')
rownames(grantham)=a(rownames(grantham))
colnames(grantham)=a(colnames(grantham))

# some code to anotate the effects of stops in dubious orfs on the essential genes  ---------------
odub=oligos[oligos$dubious, ]
ondub=oligos[!oligos$dubious,]
dubEffect = GRanges(seqnames=o.dub$chr, 
                        ranges=IRanges(start=odub$start, end=odub$end),
                        #strand=odub$codingStrand, 
                        strand=ifelse(odub$codingStrand=='+', '-', '+'),
                        unique.Index=odub$unique.Index,
                        flag.syn=odub$flagsyn,
                        id=rownames(odub), gene=odub$GENEID
                        )
replacementSeq=DNAStringSet(ifelse(odub$codingStrand=='-',  'TCA', 'TGA'))
pdub=predictCoding(dubEffect,  txdb, sacCer3, replacementSeq, ignore.strand=FALSE)

pdub.df=data.frame(pdub)

oblupdf=data.frame(oligo.blup=rfull3$oligo[,1], oligo.blup.var=attr(rfull3$oligo, 'postVar')[,1,])
rownames(oblupdf)=rownames(rfull3$oligo)

pdub.effect=merge(pdub.df,oblupdf, by.x='id', by.y='row.names')
#attr(rfull3$oligo, 'postVar')[,1,]
#weight or threshold correlations 

granth.effects=list() 
for(i in 1:length(pdub)) {
    goi=pdub$GENEID[i]
    loi=pdub$PROTEINLOC[[i]]
    varoi=s2c(as.character(pdub$VARAA[i]))
    refoi=s2c(as.character(pdub$REFAA[i]))
    idoi=pdub$id[i]
  
    ec=grep('\\*', refoi)
    if(length(ec)>0 ) {
        varoi=varoi[-ec]
        refoi=refoi[-ec]
    }
    vrlookup=cbind(refoi, varoi)
    granth.effects[[idoi]]=grantham[vrlookup]
}
}

prov.effects=list()
for(i in 1:length(pdub)) {
    goi=pdub$GENEID[i]
    loi=pdub$PROTEINLOC[[i]]
    varoi=s2c(as.character(pdub$VARAA[i]))
    idoi=pdub$id[i]
    peff=try({provean_scores[[goi]][cbind(loi,match(varoi, colnames(provean_scores[[goi]])))]})
    if(class(peff)=='try-error') { prov.effects[[idoi]]= NA } 
    else {    prov.effects[[idoi]]=provean_scores[[goi]][cbind(loi,match(varoi, colnames(provean_scores[[goi]])))] }
}

max.prov.effect=sapply(prov.effects, function(x) x[which.min((x))] )
sum.prov.effect=sapply(prov.effects, function(x) sum(x,na.rm=T) )
max.granth.effect=sapply(granth.effects, function(x) x[which.max((x))] )
sum.granth.effect=sapply(granth.effects, function(x) sum(x,na.rm=T) ) #x[which.max((x))] )

mpe=sapply(max.prov.effect, function(x){ if(length(x)==0 ){return(NA)} else{x}  })
mpe2=sapply(sum.prov.effect, function(x){ if(length(x)==0 ){return(NA)} else{x}  })
mpe3=sapply(max.granth.effect, function(x){ if(length(x)==0 ){return(NA)} else{x}  })
mpe4=sapply(sum.granth.effect, function(x){ if(length(x)==0 ){return(NA)} else{x}  })

mp.df=data.frame(id=names(mpe), max.prov=as.numeric(mpe), sum.prov=as.numeric(mpe2), max.granth=as.numeric(mpe3), sum.granth=as.numeric(mpe4), stringsAsFactors=F)
pdub.effect=merge(pdub.effect, mp.df, by='id', all.x=T)
cor.test(pdub.effect$max.prov[pdub.effect$oligo.blup.var<1.5], pdub.effect$oligo.blup[pdub.effect$oligo.blup.var<1.5], method='pearson')

ycw=cbind(pdub.effect$max.prov, pdub.effect$oligo.blup)
yw= pdub.effect$oligo.blup.var
byy=which(is.na(ycw))
ycw=ycw[-byy,]
yw=yw[-byy]


pndub.effect=merge(ondub,rfull3$oligo, by.x='row.names', by.y='row.names')

pdub.split=split(pdub.effect$'(Intercept)', pdub.effect$CONSEQUENCE)
stripchart(pdub.split, vertical=T, method='jitter', pch=21)
boxplot( pdub.split, add=T)

pdub.split2=pdub.split
pdub.split2[['essentials']]=red.effs$oligo[,1]
stripchart(pdub.split2, vertical=T, method='jitter', pch=21)
boxplot( pdub.split2, add=T)

summary(lm(scale(pdub.effect$'(Intercept)')~pdub.effect$VARAA))
saveRDS(pdub.effect, file='~/predictDubious.RDS')

pdub.no.over=(odub[!(rownames(odub) %in% pdub$id),])
pdub.no.over.effect=merge(pdub.no.over,rfull3$oligo, by.x='row.names', by.y='row.names')
saveRDS(pdub.no.over.effect, file='~/intergenic_oligos.RDS')


pdub.split[['Intergenic']]=pdub.no.over.effect$'(Intercept)'
pdub.split[['Essential']]=pndub.effect$'(Intercept)'

stripchart(pdub.split, vertical=T, method='jitter', pch=21)
boxplot( pdub.split, add=T)

x=cbind(pdub.effect$REFAA, pdub.effect$VARAA)
x1=strsplit(x[,1], '')
x2=strsplit(x[,2], '')
l1=mapply(function(x,y){x[1]!=y[1]}, x1,x2 )
l2=mapply(function(x,y){x[2]!=y[2]}, x1,x2 )

l2[is.na(l2)]=FALSE



mtest=merge(pdub.df, mp.df, by='id', all.x=T)



#pdub1=predictCoding(dubEffect,  txdb, sacCer3, reverseComplement(DNAStringSet(odub$VARCODON)), ignore.strand=FALSE)
#pdub2=predictCoding(dubEffect,  txdb, sacCer3,(DNAStringSet(odub$VARCODON)), ignore.strand=FALSE)
#
#cbind(as.character(strand(pdub1)), as.character(pdub1$VARAA),  as.character(pdub2$VARAA))
#
#choosecbind( as.character(pdub1$VARAA),  as.character(pdub2$VARAA))
#id.ref=paste(pdub$id, pdub$VARCODON, sep=':')
#o.ref=paste(rownames(odub), odub$VARCODON, sep=':')
#pdub2=pdub[!id.ref %in% o.ref]
#
#
#pds=split(pdub, pdub$id)
#pds2=list()
#for(i in names(pds)) {
#    pds2[[i]]=pds[[i]][as.character(pds[[i]]$REFAA)!=odub[i,'REFAA'],]
#}
#
#
#
#pdub=pdub[pdub$CONSEQUENCE!='nonsense',]
#
#
#          
#pdub=pdub[as.character(pdub$gene)==as.character(pdub$GENEID),]
#pdub[pdub$CONSEQUENCE=='nonsense',]

