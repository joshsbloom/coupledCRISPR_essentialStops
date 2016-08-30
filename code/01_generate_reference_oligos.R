# This script generates the expected sequences for the essential stop oligo library (post 2nd-step cloning)
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
eagI='CGGCCG'
# search for 
#more than 1 instance of 
MluI='ACGCGT'
BstEII='GGTNACC'
SphI='GCATGC' 

#original oligos location '/data/CRISPR_nonsense/essentials_choose10_101bp_final_construct_skpRF2.txt'
oligos=read.delim(file='/media/jbloom/d1/coupled_CRISPR/Reference/essentials_choose10_101bp_final_construct_skpRF2.txt', sep='\t',stringsAsFactors=F, header=T)

# beginning of read
#GCGGTATTTCACACCGCATACGTCAGATGTGTATAAGAGACAG
#NNNNN
#ACGCGT
#[Repair template]
#CCAGCTTCACA[CGGC]
#cggccggtacccaattcgccctatagtgagtcgtat
#end of read 

final.oligoDrc=DNAStringSet(oligos$final.oligo.toSynth)
# Sequence to determine if oligo needs to be reverse complemented (next time add this to annotation table for oligo before printing(!) 
revme='GCCGTGTGAAGCTGG'
flag.rc.oligos=vcountPattern(revme, final.oligoDrc)
#for(i in 1:length(final.oligoDrc) ){
#    print(i)
#    if(flag.rc.oligos[i]) final.oligoDrc[i]=reverseComplement(final.oligoDrc[i]) 
#}
final.oligoDrc[which(flag.rc.oligos==1)]=reverseComplement(final.oligoDrc[which(flag.rc.oligos==1)]) 


# for analysis with guide (after first step cloning) (this includes adapter sequences)
#chopped.oligo=DNAStringSet(as.character(final.oligoDrc), start=16,end=186)

# for analysis without guide (after second step cloning)
# goal is to extract just the repair template sequence
chopped.oligo=DNAStringSet(as.character(final.oligoDrc), start=71,end=171)
oligo.name=apply(oligos,1,function(x) paste(x[c(1:3, 5:23, 25,26)], collapse=':'))
oligo.name=gsub(' ', '', oligo.name)
names(chopped.oligo)=oligo.name
#write.fasta(chopped.oligo, file='/media/jbloom/d1/coupled_CRISPR/042916_Estops_Hiseq/eStops.fasta')
#, names=oligo.name,

# original output was here '/media/jbloom/d1/coupled_CRISPR/042916_Estops_Hiseq/ref/eStops_repair_template_only.fasta'
writeXStringSet(chopped.oligo , file='/media/jbloom/d1/coupled_CRISPR/Reference/eStops_repair_template_only.fasta')
# original output was here '/media/jbloom/d1/coupled_CRISPR/042916_Estops_Hiseq/ref/eStops_oligos.RData'
save(oligos, file='/media/jbloom/d1/coupled_CRISPR/Reference/eStops_oligos.RData')

# NOTE run this once from the command line 
#bwa index eStops_repair_template_only.fasta


#>265
#cgttcgaaacttctccgcagtgaaagata [insert here] cggtacccaattcgccctatagtgagtcgtat
#>275
#gttcgaaacttctccgcagtgaaagata [insert here] cggtacccaattcgccctatagtgagtcgtattacgcgcgctcactggccgtcgttttacaacgtcgtgactgggaaaac


