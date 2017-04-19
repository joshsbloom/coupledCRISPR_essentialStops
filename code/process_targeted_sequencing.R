# process miseq targeted sequencing of crispr edits
# modifierd 11/30/2016 to handle Olga's and Connie's multi-multiplexed libraries (multiple sites per index)

#find /media/jbloom/d1/CRISPR_base_editor/targeted_sequencing/CRISPR_Pilot2-33767777/  -type f -iname '*.fastq.gz' -exec mv {} /media/jbloom/d1/CRISPR_base_editor/targeted_sequencing/fastq/ \;
# then used thunar in linux for batch rename to parse off text between name given in csv key file and name from miseq

library(Biostrings)
library(seqinr)
library(seqLogo)
library(ShortRead)
library(rbamtools)
library("BSgenome.Scerevisiae.UCSC.sacCer3")
library(gdata)
sacCer3=BSgenome.Scerevisiae.UCSC.sacCer3

#----------------------------------------------------------------------------------------------------------------------
# change these variables as necessary------------------------------------------------------------
# put miseq data in fastq/ directory ... rename appropriately (e.g. OS_pool1_R1.fastq.gz)
git.dir='/media/jbloom/d1/coupled_CRISPR/Experiments/Preliminary/Pilot_092015/'

#Load accessory functions
source(paste0(git.dir, 'process_targeted_sequencing_accessory_functions.R'))

fastq.dir=paste0(git.dir, 'fastq/')
out.dir=paste0(git.dir, 'out/') 
dir.create(out.dir)

trimmomatic.call='java -jar /home/jbloom/Local/Trimmomatic-0.32/trimmomatic-0.32.jar SE -threads 71 ' # trailing space is necessary'
# adapters are here (check these)
adapter.file=paste0(git.dir, 'Nextera_primers_PE.fa')
# key file
#sample.key.file=paste0(git.dir, 'site_info_MS_plus.csv')
sample.key.file=paste0(git.dir, 'site_info_MS.xlsx')

bwa.call='bwa mem -t 71'
pear.params='-v 10 -m 400 -j 71 -y 8G -q 20'
#-----------------------------------------------------------------------------------------------------------


# build lookup table for samples ---------------------------------------------------------------------------
samples=buildLookupTable(sample.key.file, sacCer3)
# organize input by library
samples.by.library=split(samples, samples$library)

# ! Preprocessing steps  ---------------------------------------------------------------------------------
   
## run pear ----------------------------------------------------------------------------------------------
runPear(samples[26:27,], pear.params)

# run Trimmomatic
runTrimmomatic(samples[26:27,], trimmomatic.call, adapter.file)


# read alignment
alignReads(samples.by.library[c(9,24)], out.dir)

# read alignments into R 
#sample.alignments=parseAlignments(samples.by.library, samples, out.dir)

bamRange.list=buildBamRange_list(samples.by.library, samples, out.dir) 
sapply(bamRange.list, function(x) sapply(x,nrow) )
#save(bamRange.list, file = paste0(git.dir, 'bamRange.list.RData'))

load(paste0(git.dir, 'bamRange.list.RData'))
###############################################################################################
sample.vec=paste(samples$library, samples$samp, samples$chrom, samples$start, samples$end, sep=':')
schar=t(apply(samples[,-9],1,paste ,collapse=' '))

alignments=list()
#pdf(file=paste0(git.dir, 'mismatches2.pdf'), width=15, height=10)
threads=10
#library(foreach)
for (libr in names(samples.by.library)[c(9,24)] ){  
    print(libr)
    for(samp in names(bamRange.list[[libr]]) ) {
        print(samp)
        observed.seqs=bamRange.list[[libr]][[samp]]$seq
        fobserved.seqs=as.factor(observed.seqs)
        to.align.seqs=as.character(levels(fobserved.seqs))
        idx=match(samp, sample.vec)
        subject.seq=samples$target.fasta[idx]
            
        test= cut(1:length(to.align.seqs), threads)
        tas=split(to.align.seqs,test)
        mat = nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = FALSE, type='DNA')
        mpas=mclapply(tas, function(x) {
                     pairwiseAlignment(pattern=DNAStringSet(x), 
                                        subject=subject.seq,type='local', substitutionMatrix=mat, gapOpening=5, gapExtension=2) }, mc.cores=threads)

        ampas.mat=do.call('rbind', mclapply(mpas, function(x) as.matrix(aligned(x)), mc.cores=threads))
        tot.seqs=as.vector(table(fobserved.seqs))
        tot.seqs.vec=rep(1:length(tot.seqs),tot.seqs)
        tot.seqs.split=split(tot.seqs, test)
    
        mismatches=mcmapply(function(x,y) { 
                                m=mismatchTable(x) 
                                m$instances=y[m$PatternId]
                                return(m)
                     }
                     ,x=mpas, y=tot.seqs.split, mc.cores=threads, SIMPLIFY=FALSE)

        mismatches=do.call('rbind', mismatches)
        mismatches=mismatches[order(mismatches$instances, decreasing=T),]

        mismatch.table=sapply( split(mismatches, mismatches$SubjectStart), function(x) table(rep(x$PatternSubstring, x$instances)) )
        #crude.plot=sapply( split(mismatches, mismatches$SubjectStart), function(x) sum(x$instances)) #/length(tot.seqs.vec)
        mismatch.subject.pos=as.numeric(names(split(mismatches, mismatches$SubjectStart)))
        colnames(mismatch.table)=s2c(as.character(subject.seq))[ mismatch.subject.pos]
        attr(mismatch.table, 'chr.pos')= mismatch.subject.pos + samples[idx, 'start']
        attr(mismatch.table, 'total.reads')=length( tot.seqs.vec)

        subject.start.list=mclapply(mpas, function(x) {start(subject(x))}, mc.cores=threads)
        
        deletions=mcmapply(function(pas, subject.start, tot.seqs) {
                            deletion.intervals=RangedData(deletion(pas))
                            deletion.intervals$subject.start=subject.start[as.numeric(deletion.intervals$space)]
                            deletion.intervals$start=start(deletion.intervals)+deletion.intervals$subject.start-1
                            deletion.intervals$end= end(deletion.intervals)+deletion.intervals$subject.start-1
                            deletions=data.frame(PatternID=as.numeric(deletion.intervals$space), start=deletion.intervals$start, end=deletion.intervals$end, width=width(deletion.intervals))
                            deletions$instances=tot.seqs[deletions$PatternID]
                            deletions=deletions[order(deletions$instances, decreasing=T),]
                            return(deletions) }, pas=mpas, subject.start=subject.start.list, tot.seqs=tot.seqs.split,mc.cores=threads, SIMPLIFY=FALSE)
        deletions=do.call('rbind', deletions)
        deletions=deletions[order(deletions$instances, decreasing=T),]

        insertions=mcmapply(function(pas, subject.start, tot.seqs) {
                            insertion.intervals=RangedData(insertion(pas))
                            insertion.intervals$subject.start=subject.start[as.numeric(insertion.intervals$space)]
                            insertion.intervals$start=start(insertion.intervals)+insertion.intervals$subject.start-1
                            insertion.intervals$end= end(insertion.intervals)+insertion.intervals$subject.start-1
                            insertions=data.frame(PatternID=as.numeric(insertion.intervals$space), start=insertion.intervals$start, end=insertion.intervals$end, width=width(insertion.intervals))
                            insertions$instances=tot.seqs[insertions$PatternID]
                            insertions=insertions[order(insertions$instances, decreasing=T),]
                            return(insertions) }, pas=mpas, subject.start=subject.start.list, tot.seqs=tot.seqs.split,mc.cores=threads, SIMPLIFY=FALSE)
        insertions=do.call('rbind', insertions)
        insertions=insertions[order(insertions$instances, decreasing=T),]

        pdf(file=paste0(git.dir, make.names(samp),'.pdf'), width=15, height=10)
            par(mfrow=c(2,1), oma=c(2,2,2,2))
            plot( attr(mismatch.table, 'chr.pos'),colSums(mismatch.table), type='n', xlab=samples[idx, 'chrom'], ylab='counts')
            text(attr(mismatch.table, 'chr.pos'), mismatch.table['A',], 'A')
            text(attr(mismatch.table, 'chr.pos'), mismatch.table['C',], 'C')
            text(attr(mismatch.table, 'chr.pos'), mismatch.table['T',], 'T')
            text(attr(mismatch.table, 'chr.pos'), mismatch.table['G',], 'G')
            axis(3, at=attr(mismatch.table, 'chr.pos'), colnames(mismatch.table))
            abline(v=samples[idx,'pam'], col='red')
            plot( attr(mismatch.table, 'chr.pos'),colSums(mismatch.table), type='n', xlab=samples[idx, 'chrom'], ylab='counts', xlim=c(samples[idx,'pam']-20, samples[idx,'pam']+20))
            text(attr(mismatch.table, 'chr.pos'), mismatch.table['A',], 'A')
            text(attr(mismatch.table, 'chr.pos'), mismatch.table['C',], 'C')
            text(attr(mismatch.table, 'chr.pos'), mismatch.table['T',], 'T')
            text(attr(mismatch.table, 'chr.pos'), mismatch.table['G',], 'G')
            axis(3, at=attr(mismatch.table, 'chr.pos'), colnames(mismatch.table))
            abline(v=samples[idx,'pam'], col='red')
            title(samp,outer=T)
        dev.off()
        alignments[[libr]][[samp]]$subject.seq=subject.seq
        alignments[[libr]][[samp]]$alignment=mpas
        alignments[[libr]][[samp]]$tot.seqs=tot.seqs.split
        alignments[[libr]][[samp]]$mismatches=mismatches
        alignments[[libr]][[samp]]$deletions=deletions
        alignments[[libr]][[samp]]$insertions=insertions
        alignments[[libr]][[samp]]$mismatch.table=mismatch.table
   
    }
}
#save(alignments, file=paste0(git.dir, 'alignments.RData'))
load(paste0(git.dir, 'alignments.RData'))




haplotypes=list()
for (libr in names(alignments)){  
        print(libr)
    for(samp in names(alignments[[libr]]) ) {
        #observed.seqs=bamRange.list[[libr]][[samp]]$seq
        #fobserved.seqs=as.factor(observed.seqs)
        #to.align.seqs=as.character(levels(fobserved.seqs))
        idx=match(samp, sample.vec)
        pampos=samples[idx,'pam']

        subject.seq=alignments[[libr]][[samp]]$subject.seq
        alignment=alignments[[libr]][[samp]]$alignment
        tot.seqs.split=alignments[[libr]][[samp]]$tot.seqs
        mismatches=alignments[[libr]][[samp]]$mismatches
        deletions=alignments[[libr]][[samp]]$deletions
        insertions=alignments[[libr]][[samp]]$insertions
        mismatch.table=alignments[[libr]][[samp]]$mismatch.table
        
        total.seqs=sum(sapply(tot.seqs.split, sum))
        #print(total.seqs)

        mismatches$newID=paste0(gsub('\\.[0-9]*$', '', rownames(mismatches)), mismatches$PatternId)
        mismatches$chrStart=mismatches$SubjectStart+samples[idx, 'start']
        #ms=split(mismatches, mismatches$newID)
        #ms=ms[match(unique(mismatches$newID), names(ms))]
 
        deletions$newID=paste0(gsub('\\.[0-9]*$', '', rownames(deletions)), deletions$PatternID)
        deletions$chrStart=deletions$start+samples[idx, 'start']
        insertions$newID=paste0(gsub('\\.[0-9]*$', '', rownames(insertions)), insertions$PatternID)
        insertions$chrStart=insertions$start+samples[idx, 'start']

        dframe=data.frame(PatternId=deletions$PatternID,
                          PatternStart=0,
                          PatternEnd=0,
                          PatternSubstring='-',
                          SubjectStart=deletions$start,
                          SubjectEnd=deletions$end,
                          SubjectSubstring='-',
                          instances=deletions$instances,
                          newID=deletions$newID,
                           chrStart=deletions$chrStart, stringsAsFactors=F)
        iframe=data.frame(PatternId=insertions$PatternID,
                          PatternStart=0,
                          PatternEnd=0,
                          PatternSubstring='-',
                          SubjectStart=insertions$start,
                          SubjectEnd=insertions$end,
                          SubjectSubstring='-',
                          instances=insertions$instances,
                          newID=insertions$newID,
                           chrStart=insertions$chrStart, stringsAsFactors=F)

        m=rbind(mismatches,dframe,iframe)
        ms=split(m, as.character(m$newID))
        ms=ms[match(unique(c(mismatches$newID, dframe$newID, iframe$newID)), names(ms))]
        haplotypes[[libr]][[samp]]=ms
        }
}

#reduce haplotypes
#mshaps=sapply(ms, function(x) paste(unlist(x$PatternStart), unlist(x$PatternEnd), unlist(x$PatternSubstring), collapse=':') ) #

#save(haplotypes, file=paste0(git.dir, 'haplotypes.RData'))


load(paste0(git.dir, 'haplotypes.RData'))
#classify haplotypes
pdf(file=paste0(git.dir, 'hist_non_hap_matches_32','.pdf'), width=15, height=10)
summary.stats=list()
for (libr in names(alignments)){  
        print(libr)
    for(samp in names(alignments[[libr]]) ) {
        #observed.seqs=bamRange.list[[libr]][[samp]]$seq
        #fobserved.seqs=as.factor(observed.seqs)
        #to.align.seqs=as.character(levels(fobserved.seqs))
        idx=match(samp, sample.vec)
        pampos=samples[idx,'pam']
        ms=haplotypes[[libr]][[samp]]
        tot.seqs.split=alignments[[libr]][[samp]]$tot.seqs
        total.seqs=sum(sapply(tot.seqs.split, sum))

        # classify
        #1 = matches haplotype no mismatches or indels in 50bp window of pam
        #2 = indel in 50 bp window around pam
        #4 = WT
        #31 = partial edit (+ mismatch)
        #32 = full edit (+ mismatch)
        #33 = no edit (+ mismatch)

        #!3 = everything else
        # get index in target hap of expected hap
        expAlt=ms[[1]]$PatternSubstring
        expRef=ms[[1]]$SubjectSubstring
        expStart=ms[[1]]$SubjectStart
        
       rp=rep(NA, length(ms))
       #patterns= mclapply(ms, function(x) {
       #x=ms[[i]]
       pb =txtProgressBar(min = 1, max = length(ms), style = 3)
       for(i in 1:length(ms)) {
            setTxtProgressBar(pb, i)
            x=ms[[i]]
            cl=4
            idmatch=x$SubjectStart %in% expStart
            #x=ms[[9]]
            #full match
            other.in.pam.window=sum(x$chrStart[!idmatch]> (pampos-25)  & x$chrStart[!idmatch]< (pampos+25))
            if(other.in.pam.window>0) { cl=33 }
            del.in.pam.window=sum(x$chrStart[!idmatch & (x$PatternSubstring=='-')]> (pampos-25)  & x$chrStart[!idmatch& (x$PatternSubstring=='-')]< (pampos+25))
            if(del.in.pam.window) { cl=2 }
            match.length=sum(idmatch)
            hap.match=sum(x$PatternSubstring[idmatch]==expAlt)
            if(match.length==length(expStart)) {
                if(hap.match==length(expStart)) { cl=1 }  else {cl=31}
            }
            if(cl==1 & other.in.pam.window>0) { cl=32 }

            #partial match
            if( match.length>0  & match.length<length(expStart)) {
                if(hap.match>0) {cl=31}
            }
            rp[i]=cl
            #print(paste(i,cl))
          }
        close(pb)
        
        patterns=cbind(rp, sapply(ms, function(x)x$instances[1]))

        #rp=do.call('rbind', patterns)
        pattern.classes=sapply(split(patterns[,2],patterns[,1]), sum)
        pattern.classes['4']= pattern.classes['4']+total.seqs-sum(pattern.classes)
        #pattern.classes=c(pattern.classes, total.seqs-sum(pattern.classes))
        #names(pattern.classes)[4]='4'
        print(schar[idx])
        print(pattern.classes)


       all32=ms[which(rp==32)]
       xx2=unlist(sapply(all32, function(x) rep(x$PatternStart, x$instances)))
       xx=hist(xx2, breaks=100000,plot=F)
       xx$counts=xx$counts/sum(patterns[patterns[,1]==32,2])
       plot(xx, xlab='pos', ylab='fraction of total class 32 reads', main=samp)
       #print(pattern.classes)
       summary.stats[[samp]]=pattern.classes
    }
}
dev.off()

plyr::rbind.fill(summary.stats)
summary.counts=plyr::rbind.fill(lapply(summary.stats, function(x) data.frame(t(x))))
slabels=apply(samples[,1:8],1, paste, collapse=':')
idx=match(names(summary.stats), sample.vec)
rownames(summary.counts) = slabels[idx]
library(WriteXLS)
WriteXLS(summary.counts, '~/Desktop/summaryCountsPilot.xls', row.names=T)
summary.stats=summary.counts
summary.stats[is.na(summary.stats)]=0
save(summary.stats, file = paste0(git.dir, 'summaryStats.RData'))




#rownames(summary.counts)=names(summary.stats)
#
#load(paste0(git.dir, 'summaryStats.RData'))
#

#slabels=apply(samples[,1:8],1, paste, collapse=':')
#idx=match(names(summary.stats), sample.vec)
#names(summary.stats) = slabels[idx]
#save(summary.stats, file = paste0(git.dir, 'summaryStats.RData'))


#pdf(file=paste0(git.dir, make.names(samp),'.pdf'), width=15, height=10)
#            par(mfrow=c(2,1), oma=c(2,2,2,2))
#            plot( attr(mismatch.table, 'chr.pos'),colSums(mismatch.table), type='n', xlab=samples[idx, 'chrom'], ylab='counts')
#            text(attr(mismatch.table, 'chr.pos'), mismatch.table['A',], 'A')

#text(attr(mismatch.table, 'chr.pos'), mismatch.table['C',], 'C')
#            text(attr(mismatch.table, 'chr.pos'), mismatch.table['T',], 'T')
#            text(attr(mismatch.table, 'chr.pos'), mismatch.table['G',], 'G')
#            axis(3, at=attr(mismatch.table, 'chr.pos'), colnames(mismatch.table))
#            abline(v=samples[idx,'pam'], col='red')
#            plot( attr(mismatch.table, 'chr.pos'),colSums(mismatch.table), type='n', xlab=samples[idx, 'chrom'], ylab='counts', xlim=c(samples[idx,'pam']-20, samples[idx,'pam']+20),
#                 ylim=c(0, max(colSums(mismatch.table)[colSums(mismatch.table)<(max(colSums(mismatch.table))/4)]) ) )
#            text(attr(mismatch.table, 'chr.pos'), mismatch.table['A',], 'A')
#            text(attr(mismatch.table, 'chr.pos'), mismatch.table['C',], 'C')
#            text(attr(mismatch.table, 'chr.pos'), mismatch.table['T',], 'T')
#            text(attr(mismatch.table, 'chr.pos'), mismatch.table['G',], 'G')
#            axis(3, at=attr(mismatch.table, 'chr.pos'), colnames(mismatch.table))
#            abline(v=samples[idx,'pam'], col='red')
#            title(samp,outer=T)
 #           readline()
       # dev.off()
  return(c(cl, x$instances[1]) ) }, mc.cores=70)
             rp=do.call('rbind', patterns)



            #if(sum(idmatch)==0) {
            #    other.in.pam.window=sum(x$chrStart[!idmatch]> (pampos-25)  & x$chrStart[!idmatch]< (pampos+25))
            #    if(other.in.pam.window>0) {
            #        cl=33
            #    }
            #}
            if(sum(idmatch)>0) {
                cl=31
            }
            if(sum(idmatch)==length(expStart)){
                hapgood=all.equal(x$SubjectStart[idmatch],expStart)==T & all.equal(x$PatternSubstring[idmatch],expAlt)==T   & all.equal(x$SubjectSubstring[idmatch], expRef)==T
                if(hapgood) {  cl=1 }
            } 
            if(sum(!idmatch)>0){
                other.in.pam.window=sum(x$chrStart[!idmatch]> (pampos-25)  & x$chrStart[!idmatch]< (pampos+25))
                del.in.pam.window=sum(x$chrStart[!idmatch & (x$PatternSubstring=='-')]> (pampos-25)  & x$chrStart[!idmatch& (x$PatternSubstring=='-')]< (pampos+25))
                # this is class 32 
                if(other.in.pam.window>0) {
                    #cl=3
                    cl=32
                } else {
                    cl=33
                }
                if(del.in.pam.window>0){
                    cl=2
                }
           }
            # patterns[[i]]=(c(cl,x$instances[1]))
            return(c(cl, x$instances[1]) )# ,x$instances[1]))
         }, mc.cores=70)

