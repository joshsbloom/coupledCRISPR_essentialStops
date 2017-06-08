library(seqinr)
library(org.Sc.sgd.db)
# the easy part -------------------------------------------
lc.fasta=read.fasta('/media/jbloom/d1/coupled_CRISPR/Reference/uniprot-sCerevisiaeProteome_160404.out.fasta',forceDNAtolower = FALSE)
lc.intervals=lapply(lc.fasta, function(x) {
                   x.ind = seq(1:length(x))
                   x.lc  = grepl('[a-z]', x)
                   xrle=rle(x.lc)
                   xrle$values=as.factor(1:length(xrle$values))
                   boundary.set=split(x.ind[x.lc], inverse.rle(xrle)[x.lc])
                   t(sapply(boundary.set[sapply(boundary.set, length)>0], range))
                })
#-----------------------------------------------------------


# now the stupid names

# keep original name to double check we haven't screwed up name matching
# rename list
uniprot.IDs_from_lc.list=sapply(strsplit(names(lc.intervals), '\\|'), function(x)x[2])

#straight from here
#http://bioconductor.org/packages/release/data/annotation/manuals/org.Sc.sgd.db/man/org.Sc.sgd.db.pdf
x <- org.Sc.sgdUNIPROT
# Get the Systematic ORF IDs that are mapped to a Uniprot ID
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])

uniprot2sys=data.frame(uniprot.key=as.vector(unlist(xx)), 
           yname.key=rep(names(xx), sapply(xx, length)), 
           stringsAsFactors=F)

ygene.matches=match(uniprot.IDs_from_lc.list, uniprot2sys[,1])

for(n in 1:length(lc.intervals)) {
    attr(lc.intervals[[n]],'fasta.name' )=names(lc.intervals)[n]
    attr(lc.intervals[[n]],'uniprot' )=uniprot.IDs_from_lc.list[n]
    attr(lc.intervals[[n]],'systematic.name' )=uniprot2sys[ygene.matches[n],2] 
}

has.lowcomp=sapply(lc.intervals, ncol)>0
sapply(lc.intervals, nrow)
# now build data.frame
lc.start.stop=do.call('rbind', lc.intervals[has.lowcomp])
gene=as.character(rep(sapply(lc.intervals[has.lowcomp], function(x) attr(x,'systematic.name')), sapply(lc.intervals[has.lowcomp],nrow)))
low.complexity=data.frame(gene=gene, start=lc.start.stop[,1], stop=lc.start.stop[,2], stringsAsFactors=F)
low.complexity=na.omit(low.complexity)
attr(low.complexity, 'na.action')=NULL

save(low.complexity,file='/media/jbloom/d1/coupled_CRISPR/Reference/low.complexity.RData')
