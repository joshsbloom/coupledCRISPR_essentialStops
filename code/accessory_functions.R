library(yeastExpData)
# Build oligo data frame
# Build annotation structures 
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
flagsyn=sort(as.vector(unlist(sapply(as.character(synolis$guide), function(x)  grep(x, oligos$guide)))))
oligos$flagsyn=FALSE
oligos$flagsyn[flagsyn]=TRUE

oligos$unique.Index=1:nrow(oligos)


data(fcabundance)
fcabundance$yORF=as.character(fcabundance$yORF)
data(gfp)
gfp$yORF=as.character(gfp$yORF)
half_lives=read.delim('/media/jbloom/d1/coupled_CRISPR/Experiments/082316/ref/pnas_0605420103_SuppDataSet.txt', header=T ,sep='\t', stringsAsFactors=F)[,1:4]
names(half_lives)[1]='yORF'

abundance=merge(gfp, fcabundance, by='yORF', all.x=T)
abundance=merge(abundance, half_lives, by='yORF', all.x=T)


# some additional information
evolvability=read.delim('/media/jbloom/d1/coupled_CRISPR/Reference/Evolvability.csv', header=T, sep='\t', stringsAsFactors=F)
evolvability$Systematic.name 
evolvability$X
comp_human=read.delim('/media/jbloom/d1/coupled_CRISPR/Reference/Marcotte_complementsHuman.csv', header=T, sep='\t', stringsAsFactors=F)
comp_human$ScENSP
comp_human$Final.CompStatus


#haploinsufficiency
haploinf=read.delim('/media/jbloom/d1/coupled_CRISPR/Reference/haploinsufficiency.txt',  header=T ,sep='\t', stringsAsFactors=F)
haploinf$orf=gsub(' ','',haploinf$orf)
haploinf$HET_AV=as.numeric(haploinf$HET_AV)
haploinf$HET_AV=as.numeric(haploinf$HET_AV)

# set up for GO enrichment analysis -----------------------------------------------------------------------------------------------
#gcoord.key= build.gcoord.key('/data/eQTL/reference/sacCer3.fasta')

# setup for GO enrichment analysis ------------------------------------------------------------------------
# slightly modified by JB, original code from Frank Albert
myGene2GO.full.table=read.delim('/data/eQTL/reference/go_annotations_sgd.txt',  sep='\t', header=F, stringsAsFactors=F)

# data frame of gene then GO
myGene2GO=cbind( sapply(strsplit(myGene2GO.full.table[,11], '\\|'), function(x) x[1]), myGene2GO.full.table[,5])
myGene2GO=na.omit(myGene2GO)

SYS2ORF=unique(cbind(sapply(strsplit(myGene2GO.full.table[,11], '\\|'), function(x) x[1]), myGene2GO.full.table[,3]))
SYS2ORF.key=list()
SYS2ORF.key[SYS2ORF[,1]]=SYS2ORF[,2]
gene2GOList = lapply(unique(myGene2GO[,1]), function(x){myGene2GO[myGene2GO[,1] == x, 2]})
names(gene2GOList) = unique(myGene2GO[,1])

GO2geneList = lapply(unique(myGene2GO[,2]), function(x){myGene2GO[myGene2GO[,2] == x, 1]})
names(GO2geneList) = unique(myGene2GO[,2])


plotGOToTree <- function(GOdat, GOres, sigThres = 0.0005){
     # only plot if there are any significant GO terms (SEE ABOVE for 
    #"significance"; I am somewhat lenient here):
    #     # we need these extra lines because very small p-values are 
    #reported as a text string "< X", rather than a numeric
     toTest <- as.numeric(GenTable(GOdat, pVal = GOres)[1,6])
     if(is.na(toTest)){toTest <- 0.000000000000000000000000000001}
     if (toTest < sigThres){
         showSigOfNodes(GOdat, score(GOres), firstSigNodes = 
        length(score(GOres)[score(GOres) < sigThres]), useInfo = "all")
     }else{
         plot(1,1)
         legend("topright", legend="no significant GO categories", 
        cex=0.4, box.lty=0)
     }
}

