library("BSgenome.Scerevisiae.UCSC.sacCer3")
sacCer3=BSgenome.Scerevisiae.UCSC.sacCer3
library(seqinr)
library(intervals)
library("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene")
txdb=TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
library(VariantAnnotation)
library(foreach)
library(doMC)
registerDoMC(cores=6)

PAMp_string=DNAString('GG')
PAMm_string=reverseComplement(PAMp_string)

registerDoMC(cores=5)

exonsByGene=(exonsBy(txdb, by='gene'))
protein.lengths=lapply(width(exonsByGene), sum)


skpF1='GGGTCACGCGTAGGA' #3’  sharon2012_F
skpF2='CGCGTCGAGTAGGGT' #3’  firstAUG_F
skpF3='CGATCGCCCTTGGTG' #3’  firstUTR_F
skpF4='GGTCGAGCCGGAACT' #3’  upstream_F
skpF5='TCCCGGCGTTGTCCT'
skpF6='CGCAGGGTCCAGAGT'

skpR1r='GTTCCGCAGCCACAC' #3’ sharon2012_R
skpR2r='GCCGTGTGAAGCTGG' #3’ firstAUG_R
skpR3r='GGTTTAGCCGGCGTG' #3’ firstUTR_R
skpR4r='GGATGCGCACCCAGA' #3’ upstream_R
skpR5r='GCTCCGTCACTGCCC'
skpR6r='GTTCGCGCGAAGGAA'

skpR1=reverseComplement(DNAString(skpR1r))
skpR2=reverseComplement(DNAString(skpR2r))
skpR3=reverseComplement(DNAString(skpR3r))
skpR4=reverseComplement(DNAString(skpR4r))
skpR5=reverseComplement(DNAString(skpR5r))
skpR6=reverseComplement(DNAString(skpR6r))

#search for 
#1 or more instance of 
eagI='CGGCCG'
# search for 
#more than 1 instance of 
MluI='ACGCGT'
BstEII='GGTNACC'
SphI='GCATGC' 



build.oligos=function(chr,txdb, PAMp, PAMm, target.genes )  { 
        print(chr)
        # Build genomic ranges of PAM sites
        gPAM = GRanges(seqnames=names(sacCer3)[chr], 
                       ranges=IRanges(start=c(start(ranges(PAMp))-1, start(ranges(PAMm))), 
                                           end=c(end(ranges(PAMp)), end(ranges(PAMm))+1)),
        strand=c(rep('+', length(PAMp)),rep('-', length(PAMm)) ) )
        genome(gPAM)='sacCer3'
        
        # Check if PAM site is coding
        codingPAM=locateVariants(gPAM, txdb, CodingVariants(), ignore.strand=TRUE)
        
        # pull out target genes 
        all.e.PAM=codingPAM$GENEID %in% target.genes
        all.coding.PAM=c(ranges(codingPAM[all.e.PAM]))

        codingpams=match(ranges(gPAM), all.coding.PAM)
        gPAM$coding=ifelse(is.na(codingpams),FALSE,TRUE)
        gPAM$Index=1:length(gPAM)
        gPAM=gPAM[gPAM$coding==1]

        gPAM$guide=rep('', length(gPAM))
        startg=as.numeric(start(gPAM))
        endg=as.numeric(end(gPAM))
        strands=as.character(strand(gPAM))
       
        gPAM_p=gPAM
        start(gPAM_p)=start(gPAM)-20
        # +2 if include PAM
        end(gPAM_p)=start(gPAM)-1
        gPAM_p=gPAM_p[strands=='+']
        gPAM$guide[strands=='+']=as.vector(as.character(getSeq(sacCer3, gPAM_p)))

        gPAM_m=gPAM
        #  -2 if include PAM 
        start(gPAM_m)=end(gPAM)+1
        end(gPAM_m)=end(gPAM)+20
        gPAM_m=gPAM_m[strands=='-']
        gPAM$guide[strands=='-']=as.vector(as.character((getSeq(sacCer3, gPAM_m))))

        nPAMs=list()
        ii=1
        for(i in c(0,1,2) ) {
            gPAM.u1=gPAM[which(strand(gPAM)=='+')]
            gPAM.u1$guideStrand=rep("+", length(gPAM.u1))
            gPAM.u1$case=rep(i, length(gPAM.u1))
            start(gPAM.u1)=start(gPAM.u1)+i
            end(gPAM.u1)=end(gPAM.u1)+i

            strand(gPAM.u1)=rep("+", length(gPAM.u1))
            gPAM.u1$codingStrand=rep("+", length(gPAM.u1))
            pc.1=predictCoding(gPAM.u1, txdb, sacCer3, DNAStringSet(rep('TGA', length(gPAM.u1))), ignore.strand=FALSE)
            nPAMs[[ii]]=pc.1[pc.1$CONSEQUENCE=='nonsense']
            ii=ii+1
            
            strand(gPAM.u1)=rep("-", length(gPAM.u1))
            gPAM.u1$codingStrand=rep("-", length(gPAM.u1))
            pc.1=predictCoding(gPAM.u1, txdb, sacCer3, DNAStringSet(rep('TCA', length(gPAM.u1))), ignore.strand=FALSE)
            nPAMs[[ii]]=pc.1[pc.1$CONSEQUENCE=='nonsense']
            ii=ii+1
        }
        for(i in c(0,-1,-2) ) {
            gPAM.u1=gPAM[which(strand(gPAM)=='-')]
            gPAM.u1$guideStrand=rep("-", length(gPAM.u1))
            gPAM.u1$case=rep(i, length(gPAM.u1))
            start(gPAM.u1)=start(gPAM.u1)+i
            end(gPAM.u1)=end(gPAM.u1)+i

            strand(gPAM.u1)=rep("+", length(gPAM.u1))
            gPAM.u1$codingStrand=rep("+", length(gPAM.u1))
            pc.1=predictCoding(gPAM.u1, txdb, sacCer3, DNAStringSet(rep('TGA', length(gPAM.u1))), ignore.strand=FALSE)
            nPAMs[[ii]]=pc.1[pc.1$CONSEQUENCE=='nonsense']
            ii=ii+1
            
            strand(gPAM.u1)=rep("-", length(gPAM.u1))
            gPAM.u1$codingStrand=rep("-", length(gPAM.u1))
            pc.1=predictCoding(gPAM.u1, txdb, sacCer3, DNAStringSet(rep('TCA', length(gPAM.u1))), ignore.strand=FALSE)
            nPAMs[[ii]]=pc.1[pc.1$CONSEQUENCE=='nonsense']
            ii=ii+1
        }

    nPAMs=do.call('c', nPAMs)
    nPAMs=nPAMs[order(start(nPAMs)),]


    nPAMs.flank=nPAMs
    #for 100 bp _oligos
    start(nPAMs.flank)=start(nPAMs)-49
    end(nPAMs.flank)=start(nPAMs)-1
    strand(nPAMs.flank)=rep('+', length(nPAMs.flank))
    refseqL=getSeq(sacCer3,nPAMs.flank)

    refseqStop=ifelse(as.vector(nPAMs$codingStrand)=='+', 'TGA', 'TCA')

    nPAMs.flank=nPAMs
    start(nPAMs.flank)=end(nPAMs)+1
    end(nPAMs.flank)=end(nPAMs)+49
    strand(nPAMs.flank)=rep('+', length(nPAMs.flank))
    refseqR=getSeq(sacCer3,nPAMs.flank)

    stopPAM_oligos=paste(refseqL, refseqStop, refseqR, sep='')
    nPAMs$stopPAM_oligos=stopPAM_oligos

    nPAMs$CDS_length=as.vector(unlist(protein.lengths[nPAMs$GENEID]))/3
    nPAMs$dist_from_CDS_end=nPAMs$CDS_length-unlist(nPAMs$PROTEINLOC)

    un=unique(nPAMs)
    df=DataFrame(un)
    dfn=data.frame(chr=as.vector(seqnames(un)),start=start(un), end=end(un),df[,-1])
    return(dfn)
}


#http://www-sequence.stanford.edu/group/yeast_deletion_project/Essential_ORFs.txt
e_orfs=(read.delim('/data/CRISPR_nonsense/E_orfs_column.txt', stringsAsFactors=F, header=F, sep='\t'))[,1]
eStops=foreach(chr=1:16, .combine=rbind) %do% {
    # Build genomic ranges of PAM sites
    PAMp=matchPattern(PAMp_string, sacCer3[[chr]], fixed=T)
    PAMm=matchPattern(PAMm_string, sacCer3[[chr]], fixed=T)
    #eStops[[chr]]=
    return(build.oligos(chr,txdb, PAMp, PAMm, e_orfs))
}
write.table(eStops, file='/data/CRISPR_nonsense/essential_stop_oligos_BY_101bp.txt',  sep='\t', quote=F, row.names=F)
save(eStops, file='/data/CRISPR_nonsense/estops_101.RData')

# Code to pick 10 oligos per essential gene weighted by distance from end of coding sequence  --------------------------------
e.by.G=split(eStops, eStops$GENEID)
# choose 10 oligos per essential gene 
e_choose10=list()
for(g in names(e.by.G) ){
        y=e.by.G[[g]]  
        x=y$dist_from_CDS_end
        l.tot=y$CDS_length
        x2=l.tot-x
        if(length(x2)>10){
            grab.inds=sort(sample(x,10, prob=x2/l.tot))
            e_choose10[[g]]=y[match(grab.inds, x),]
        }
      }
edf=do.call('rbind', e_choose10)

#----------------------------------------------------------------------------------------------------------------------------------
final.oligo=paste0(skpF2, 'GGTGACC', edf$guide,'GTTTTAGAGCATGC', 'CGATCGAT','ACGCGT', edf$stopPAM_oligos,  skpR2)
final.oligoD= DNAStringSet(final.oligo)
Acnt=letterFrequency(final.oligoD, 'A')
Tcnt=letterFrequency(final.oligoD, 'T')

final.oligoDrc=final.oligoD
for(i in 1:length(final.oligoD) ){
    print(i)
    if(Acnt[i]>Tcnt[i]) final.oligoDrc[i]=reverseComplement(final.oligoDrc[i]) 
}

eagI.cnt=vcountPattern(eagI, final.oligoDrc)
MluI.cnt=vcountPattern(MluI, final.oligoDrc)
SphI.cnt=vcountPattern(SphI, final.oligoDrc)
BstEII.cnt=vcountPattern(BstEII, final.oligoDrc, fixed=F)
bad.oligos=(eagI.cnt>0 | MluI.cnt>1 | SphI.cnt >1 | BstEII.cnt >1)

final.oligo.toSynth=final.oligoDrc
fos=data.frame(edf, final.oligo.toSynth=as.character(final.oligoDrc), stringsAsFactors=F)
fos=fos[-which(bad.oligos),]
write.table(fos,  file='/data/CRISPR_nonsense/essentials_choose10_101bp_final_construct_skpRF2.txt', sep='\t', quote=F, row.names=F)





#-----------------------------------------------
ne_orfs=unique(names(exonsByGene)[!names(exonsByGene) %in% e_orfs])
neStops=foreach(chr=1:16, .combine=rbind) %do% {
    # Build genomic ranges of PAM sites
    PAMp=matchPattern(PAMp_string, sacCer3[[chr]], fixed=T)
    PAMm=matchPattern(PAMm_string, sacCer3[[chr]], fixed=T)
    #eStops[[chr]]=
    return(build.oligos(chr,txdb, PAMp, PAMm, ne_orfs))
}
write.table(neStops, file='/data/CRISPR_nonsense/non_essential_stop_oligos_BY_101bp.txt',  sep='\t', quote=F, row.names=F)
save(neStops, file='/data/CRISPR_nonsense/nestops_101.RData')


ne.by.G=split(neStops, neStops$GENEID)
# choose 10 oligos per essential gene 
ne_choose10=list()
for(g in names(ne.by.G) ){
        y=ne.by.G[[g]]  
        grab.inds= which.min(abs(50-as.vector(unlist(as.vector(y$PROTEINLOC)))))
        #findInterval(50, as.vector(unlist(as.vector(y$PROTEINLOC))))
        ne_choose10[[g]]=y[grab.inds,]
      }
nedf=do.call('rbind', ne_choose10)

#----------------------------------------------------------------------------------------------------------------------------------
nfinal.oligo=paste0(skpF3, 'GGTGACC', nedf$guide,'GTTTTAGAGCATGC', 'CGATCGAT','ACGCGT', nedf$stopPAM_oligos,  skpR3)
nfinal.oligoD= DNAStringSet(nfinal.oligo)
Acnt=letterFrequency(nfinal.oligoD, 'A')
Tcnt=letterFrequency(nfinal.oligoD, 'T')

nfinal.oligoDrc=nfinal.oligoD
for(i in 1:length(nfinal.oligoD) ){
    print(i)
    if(Acnt[i]>Tcnt[i]) nfinal.oligoDrc[i]=reverseComplement(nfinal.oligoDrc[i]) 
}

eagI.cnt=vcountPattern(eagI, nfinal.oligoDrc)
MluI.cnt=vcountPattern(MluI, nfinal.oligoDrc)
SphI.cnt=vcountPattern(SphI, nfinal.oligoDrc)
BstEII.cnt=vcountPattern(BstEII, nfinal.oligoDrc, fixed=F)
bad.oligos=(eagI.cnt>0 | MluI.cnt>1 | SphI.cnt >1 | BstEII.cnt >1)

nfinal.oligo.toSynth=nfinal.oligoDrc
nfos=data.frame(nedf, final.oligo.toSynth=as.character(nfinal.oligoDrc), stringsAsFactors=F)
nfos=nfos[-which(bad.oligos),]
write.table(nfos,  file='/data/CRISPR_nonsense/non_essentials_choose1_50aa_101bp_final_construct_skpRF3.txt', sep='\t', quote=F, row.names=F)



ne.by.G2=split(neStops, neStops$GENEID)
# choose 10 oligos per essential gene 
ne_choose102=list()
for(g in names(ne.by.G2) ){
        y=ne.by.G2[[g]]  
        dist.vect=abs(50-as.vector(unlist(as.vector(y$PROTEINLOC))))
        first.min=which.min(dist.vect)
        dist.vect[first.min]=1e5
        grab.inds=which.min(dist.vect)
        #findInterval(50, as.vector(unlist(as.vector(y$PROTEINLOC))))
        ne_choose102[[g]]=y[grab.inds,]
      }
nedf2=do.call('rbind', ne_choose102)

#----------------------------------------------------------------------------------------------------------------------------------
nfinal.oligo2=paste0(skpF4, 'GGTGACC', nedf2$guide,'GTTTTAGAGCATGC', 'CGATCGAT','ACGCGT', nedf2$stopPAM_oligos,  skpR4)
nfinal.oligoD2= DNAStringSet(nfinal.oligo2)
Acnt=letterFrequency(nfinal.oligoD2, 'A')
Tcnt=letterFrequency(nfinal.oligoD2, 'T')

nfinal.oligoDrc2=nfinal.oligoD2
for(i in 1:length(nfinal.oligoD2) ){
    print(i)
    if(Acnt[i]>Tcnt[i]) nfinal.oligoDrc2[i]=reverseComplement(nfinal.oligoDrc2[i]) 
}

eagI.cnt=vcountPattern(eagI, nfinal.oligoDrc2)
MluI.cnt=vcountPattern(MluI, nfinal.oligoDrc2)
SphI.cnt=vcountPattern(SphI, nfinal.oligoDrc2)
BstEII.cnt=vcountPattern(BstEII, nfinal.oligoDrc2, fixed=F)
bad.oligos=(eagI.cnt>0 | MluI.cnt>1 | SphI.cnt >1 | BstEII.cnt >1)

nfinal.oligo.toSynth=nfinal.oligoDrc2
nfos2=data.frame(nedf2, final.oligo.toSynth=as.character(nfinal.oligoDrc2), stringsAsFactors=F)
nfos2=nfos2[-which(bad.oligos),]
write.table(nfos2,  file='/data/CRISPR_nonsense/non_essentials_choose2_50aa_101bp_final_construct_skpRF4.txt', sep='\t', quote=F, row.names=F)






#sampe( e.by.G[[1]]$CDS_length[1] )

#Sample 10 oligos from each essential gene starting from 



#meru_nf_orfs=c('YPL187W', 'YKL178C', 'YHR014W', 'YGL183C', 'YHL022C', 'YLL017W', 'YDL227C', 'YFR025C', 'YIL116W', 'YGL167C' )
#mStops=foreach(chr=1:16, .combine=rbind) %do% {
#    PAMp=matchPattern(PAMp_string, sacCer3[[chr]], fixed=T)
#    PAMm=matchPattern(PAMm_string, sacCer3[[chr]], fixed=T)
#    #eStops[[chr]]=
#    return(build.oligos(chr,txdb, PAMp, PAMm, meru_nf_orfs))
#}
#write.table(mStops,  file='/data/CRISPR_nonsense/Meru_target_NE_stop_oligos_BY_100bp.txt', sep='\t', quote=F, row.names=F)

