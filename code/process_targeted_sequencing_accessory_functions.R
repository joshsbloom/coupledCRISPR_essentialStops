# helper functions -----------------------------------------------------------------------------------------------------
convertPear2FASTA=function(pear.file, pear.fasta.file) {
    pearIn=scan(pear.file, what=character())
    pearSeqNames=pearIn[seq(1,length(pearIn),2)]
    pearSeqNames=paste(seq(1, length(pearSeqNames)), pearSeqNames, sep=':')
    pearSeqs=DNAStringSet(pearIn[seq(2,length(pearIn),2)])
    names(pearSeqs)=pearSeqNames
    writeXStringSet(pearSeqs,pear.fasta.file)
}

runBlat=function(blat.path='/home/jbloom/Local/bin/blat', target.fasta.file, pear.fasta.file, blat.table.file) {
    system(paste(blat.path, target.fasta.file, pear.fasta.file, blat.table.file))
}
readBlat=function(blat.table.file, nCount.filter=5){
    blat.table=read.delim(blat.table.file, header=F, sep='\t', skip=5)
    names(blat.table)=c("matches","misMatches","repMatches","nCount","qNumInsert",
                        "qBaseInsert","tNumInsert","tBaseInsert","strand","qName","qSize","qStart","qEnd",
                        "tName","tSize","tStart","tEnd","blockCount","blockSizes","qStarts","tStarts")
    blat.table$seq.counts=as.numeric(do.call('rbind', strsplit(as.character(blat.table$qName), ':'))[,2])
    blat.table=blat.table[blat.table$nCount<nCount.filter,]
    return(blat.table)
}
makeSeqLogo=function(snpConsensus) {
    #a1=floor(snpConsensus['N',]/4)
    #t1=a1+(snpConsensus['N',]-a1*4)
    snpConsensus2=snpConsensus[-5,]
    snpConsensus2[1,]=snpConsensus2[1,] #+a1
    snpConsensus2[2,]=snpConsensus2[2,] #+a1
    snpConsensus2[3,]=snpConsensus2[3,] #+a1
    snpConsensus2[4,]=snpConsensus2[4,] #+t1
    snpConsensus2=t(t(snpConsensus2)/(colSums(snpConsensus2)))
    pwm2=makePWM((snpConsensus2))
    # regions of interest
   # ic=pwm2@ic
  #  ic[1:50]=NA
  #  ic[(nchar(target.fasta)-50):(nchar(target.fasta))]=NA
  #  max.un=which.min(ic)
  #  print(max.un)
    #seqLogo(pwm2, ic.scale=F)
   # pwm3=makePWM(snpConsensus2[,(max.un-16):(max.un+16)])
    seqLogo(pwm2, ic.scale=F, xaxis=T)
    return(pwm2)
}

# build lookup table
buildLookupTable=function(sample.key.file, sacCer3) {
    #samples=read.delim(sample.key.file, header=T, sep='\t', stringsAsFactors=F)
    samples=read.xls(sample.key.file, header=T,  stringsAsFactors=F)
    chroms =paste0('chr', as.roman(seq(1:16)))
    samples$chrom=chroms[samples$chrom]
    samples= DataFrame(samples)
    #targets=list()
    samples$target.fasta=DNAStringSet(apply(samples,1, function(i) {DNAString(sacCer3[[i[3]]][i['start']:i['end']]) }))
    names(samples$target.fasta)=paste(samples$library, samples$samp, samples$chrom, samples$start, samples$end, sep=':')
    return(samples)
}

# run pear read assembler
runPear=function(samples, pear.params) {
    for(n in unique(samples$library) ) {
        # run pear
        ffq=paste0(fastq.dir, n, '_R1.fastq.gz')  #_t_R1.fastq.gz
        rfq=paste0(fastq.dir, n, '_R2.fastq.gz')  #_t_R1.fastq.gz
        ofq=paste0(fastq.dir, n, '_M') #${fastq_dir}${samp}_M
        # modfiy as appropriate
        system(paste('/home/jbloom/Local/bin/pear -f', ffq, '-r', rfq, '-o', ofq  ,pear.params))
    }
}  

# run trimmomatic 
runTrimmomatic=function(samples, trimmomatic.call, adapter.file) { 
    for(n in unique(samples$library) ) {
        m.assembled.filt =paste0(fastq.dir, n, '_M.assembled.fastq')
        mt.assembled.filt=paste0(fastq.dir, n, '_M_T.assembled.fastq')
        # run in single end mode becase reads have already been assembled with pear
        # change tri
        call.trimmomatic=paste0(trimmomatic.call,
              m.assembled.filt, ' ',
              mt.assembled.filt, ' ',
              'ILLUMINACLIP:',adapter.file,':2:20:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:200 TOPHRED33')
        system(call.trimmomatic)
    }
} 
#how to run cudadapt (need to modify to have proper adapter sequences)
#seqL1='CCTTCTTGACGAGTTCTTCTGAATTATTAACGCTTACAATTTCCTGATGCGGTATTTTCTCCTTACGCATCTGTGCGGTATTTCACACCGCTCGGAGGGCAACTATCC'
#seqR1='CCAGCTTCACACGGC'
#foreach(rr = readID)  %dopar% {
#    fin=paste(out.base.dir, 'pear/', rr, '_pear_out.assembled.fastq', sep='')
#    fout=paste(out.base.dir, 'pear/', rr, '_adapter_trimmed.fastq', sep='')
#    system(paste('~/.local/bin/cutadapt  -n 3 -g', seqL1, '-b', seqR1, '-m 40 -o', fout, fin))
#}


# run bwa and sam2pwairwise and output alignments
alignReads = function(samples.by.library, out.dir) {
    for(libr in names(samples.by.library) ){
  
        # build and index the set of expected reference sequences ------------------------
        target.fasta.file=paste0(out.dir, libr, '_ref.fasta')
        writeXStringSet(samples.by.library[[libr]]$target.fasta, target.fasta.file)
        system(paste('bwa index', target.fasta.file))
        #---------------------------------------------------------------------------------

        # align with bwa 
        mt.assembled.filt =paste0(fastq.dir, libr, '_M_T.assembled.fastq')
       sam.out.file=paste0(out.dir, libr, '_out')
      
       # generate bam alignment files
        system(paste(bwa.call, target.fasta.file, mt.assembled.filt , ' | samtools view -hSu - | samtools sort - ', sam.out.file))
      
       # specify output file location
        pretty.alignment.file=paste0(out.dir, libr, '_align.txt.gz')
        bam.file=paste0(sam.out.file, '.bam')
      
        system(paste('samtools index', bam.file)) #,  '| sam2pairwise | gzip > ', pretty.alignment.file))
         
        # filter out the reads that didn't not match any expected sequence  (q1) and convert to multiple alignment
        ## commented out this line
        #system(paste('samtools view -q1', bam.file,  '| sam2pairwise | gzip > ', pretty.alignment.file))
        
        # syntax to find and rank unmapped reads (run from shell)
        # samtools view -f4 OS_pool2_out.bam | cut -f 10 | sort -k 1,1 | uniq -c | sort -k 1,1nr > OS2_unmapped.txt
    }
}


buildBamRange_list = function(samples.by.library, samples, out.dir) {
 
    bamRange.list=list()
    for (libr in names(samples.by.library) ) {
        #for(samp in  samples.by.library[[libr]]$samp) {
         #libr='106-1'
         print(libr)
         # note that if there are more samples 
         bam.file=paste(git.dir, '/out/', libr, '_out.bam', sep='')

         # use rbamtools to extract aligned results
         reader = bamReader(bam.file,idx=TRUE)
         contigs = (refSeqDict(getHeaderText(reader)))@SN
         contigs.length=(refSeqDict(getHeaderText(reader)))@LN

          # get counts for each expected contig
          ac=bamCountAll(reader)
          bamClose(reader)

          # read in again and extract sequence data
          reader=bamReader(bam.file,idx=TRUE)
            branges=list()
            #pb =txtProgressBar(min = 1, max =length(contigs), style = 3)
            for(n in 1:length(contigs)) {
                #setTxtProgressBar(pb, n)
                # argh, 0 indexing
                branges[[contigs[n]]]=bamRange(reader, coords=c(n-1, 1, contigs.length[n]))
            }
            #close(pb)
            bamClose(reader)

            # convert the bamRange object into a list of data frames

            bdf=mclapply(branges, function(x) {
                   a=as.data.frame(x)  
                   # modified 062916
                   # modified 082416 for hiseq
                   #a$barcode=sapply(strsplit(a$name, ':'), function(x)x[2])
                   #a$barcode.type=sapply(strsplit(a$name, ':'), function(x)x[3])
                   return(a)
            }, mc.cores=4)
            # custom for one sample per library
            #if(length(bdf)==1) {   bamRange.list[[libr]]=bdf[[1]] }
            #else {
            bamRange.list[[libr]]=bdf
            #}
        #}
    }

    return(bamRange.list)

}



parseAlignments=function(samples.by.library, samples, out.dir) {
    
    alignments=list()
    for(libr in names(samples.by.library) ){
        # read in alignments
        pretty.alignment.file=paste0(out.dir, libr, '_align.txt.gz')
        fileCon = gzfile(pretty.alignment.file, open = "r")
        # note, do not need to read in ref sequence (specify as NULL)
        afilein = scan(fileCon, comment.char='', multi.line=T,
                     what=list(name='character', 
                               query='character',
                               align.stat='character',
                               ref='character'), sep="\n")
        close(fileCon)
        alignments[[libr]]=afilein
    }
    
    sample.alignments=list() 
    for(libr in names(samples.by.library) ){
        print(libr)
        afilein=alignments[[libr]]
        explode.name=strsplit(afilein$name, '\t')
        seq.name=sapply(explode.name, function(x) x[3])

        expected.samples=paste(samples$library, samples$samp, samples$chrom, samples$start, samples$end, sep=':')
        seq.match=match(seq.name, expected.samples)
        seq.match=expected.samples[seq.match]
        
        alignments.by.expected=split(afilein$query, seq.match)
         alignstat.by.expected=split(afilein$align.stat, seq.match)
          alignref.by.expected=split(afilein$ref, seq.match)
       
        for(samp in names(alignments.by.expected) ) {
                print(samp)
                que=alignments.by.expected[[samp]]
                print(length(que))
                ref=alignref.by.expected[[samp]]
                sta=alignstat.by.expected[[samp]]
                
                oq=order(as.factor(que))
                que=que[oq]
                ref=ref[oq]
                sta=sta[oq]
                
                qsr=paste(que,sta,ref, sep=":")            

                rr=rle(qsr)
                rl=rr$lengths

                rr$lengths=rr$lengths[order(rl,decreasing=T)]
                rr$values=rr$values[order(rl,decreasing=T)]

                rds=do.call('rbind', strsplit(rr$values, ':'))
                colnames(rds)=c('query', 'stat', 'ref')
                sample.alignments[[libr]][[samp]]=data.frame(rr$lengths, rds, stringsAsFactors=F)
        
        }
    }
    return(sample.alignments)
}

#make summary pdf (shows counts and  cummulative counts of most abundant sequences)
makeSummaryPDF=function(out.dir, sample.alignments) {
    pdf(file=paste0(out.dir, 'summary_stats.pdf') ,width=8, height=10)
    for(libr in names(sample.alignments) ){
        for(samp in names(sample.alignments[[libr]]) ){
             print(libr)
             print(samp)
        o2.sample=sample.alignments[[libr]][[samp]]
        par(mfrow=c(2,1))
        barplot(o2.sample[,1], names=seq(1:nrow(o2.sample)),xlim=c(1,30),
                main=paste(libr, samp), ylab='counts', sub=sum(o2.sample[,1])
                )
        plot(cumsum(o2.sample[,1])/sum(o2.sample[,1]),
         xlim=c(1,30),
         ylim=c(0,1),
         type='b',
         ylab='cummulative fraction')
        }
    }
    dev.off()
}

# output alignments as text files with frequency
makeAlignmentTxtOut=function(out.dir, sample.alignments) {
    for(libr in names(sample.alignments) ){
        for(samp in names(sample.alignments[[libr]]) ){
             print(libr)
             print(samp)
             o2.sample=sample.alignments[[libr]][[samp]]
           
             output.alignment.file=paste0(out.dir, libr, ' ', samp,  '_processed_alignment.txt')

            fileCon = file(output.alignment.file, open = "wt")
            writeLines(paste0(sprintf('%8d',o2.sample$rr.lengths), '\t', o2.sample$query, '\n',
                              sprintf('%8s', ''), '\t', o2.sample$stat, '\n',
                              sprintf('%8s', ''), '\t', o2.sample$ref, '\n'), con=fileCon)
          close(fileCon)
        }
    }
}

# makes sequence logos and returns consensus matrices
makeSequenceLogo=function(out.dir, samples, sample.alignments) {
    pdf(file=paste0(out.dir, 'sequence_logos.pdf') ,width=50, height=5)
    expected.samples=paste(samples$library, samples$samp, samples$chrom, samples$start, samples$end, sep=':')
    consensusMatrices=list()
    for(libr in names(sample.alignments) ){
        for(samp in names(sample.alignments[[libr]]) ){
        
            o2.sample=sample.alignments[[libr]][[samp]]
             print(libr)
             print(samp)

        s.ind=match(samp, expected.samples)
       
        no.clip.ref.seqs = o2.sample$ref %in% samples$target.fasta[s.ind]
        if(sum( no.clip.ref.seqs) >0 ) {
            tot.no.clip=sum(o2.sample[no.clip.ref.seqs,1])
            tot.seq=sum(o2.sample[,1])
            print(tot.no.clip)
            print(tot.seq)
            print(tot.no.clip/tot.seq)

            nc.q=o2.sample$query[no.clip.ref.seqs]
            ncq.set=DNAStringSet(nc.q)

            stacked_seq=rep(ncq.set, o2.sample[no.clip.ref.seqs,1])

            # count up bases at each position
            snpConsensus=consensusMatrix(stacked_seq, baseOnly=TRUE)
            rownames(snpConsensus)[5]='N'
            colnames(snpConsensus)=s2c(as.character(samples$target.fasta[s.ind]))
            consensusMatrices[[libr]][[samp]]=snpConsensus
            makeSeqLogo(snpConsensus) 
            }
        }
    }
    dev.off()
    return(consensusMatrices)
}



#extracting every 2nd line from 4 line chunks
#    awk 'NR%4==2' $inf > $tmp1







# blat output psl format:
#http://uswest.ensembl.org/info/website/upload/psl.html?redirect=no
#    matches    -    Number of matching bases that aren't repeats.
#    misMatches - Number of bases that don't match.
#    repMatches - Number of matching bases that are part of repeats.
#    nCount     - Number of 'N' bases.
#    qNumInsert - Number of inserts in query.
#    qBaseInsert- Number of bases inserted into query.
#    tNumInsert - Number of inserts in target.
#    tBaseInsert- Number of bases inserted into target.
#    strand     - defined as + (forward) or - (reverse) for query strand. In mouse, a second '+' or '-' indecates genomic strand.
#    qName      - Query sequence name.
#    qSize      - Query sequence size.
#    qStart     - Alignment start position in query.
#    qEnd       - Alignment end position in query.
#    tName      - Target sequence name.
#    tSize      - Target sequence size.
#    tStart     - Alignment start position in query.
#    tEnd       - Alignment end position in query.
#    blockCount - Number of blocks in the alignment.
#    blockSizes - Comma-separated list of sizes of each block.
#    qStarts    - Comma-separated list of start position of each block in query.
#    tStarts    - Comma-separated list of start position of each block in target.
#/blat 17T_S6_target.FA fastaTest.fa -out=psl out.pslu








#pdf(file=paste0(base.dir, 'plot.pdf'), width=15, height=8)
#   #targ=names(targets)[1] 
#    #target.fasta=targets$target.fasta
#    #target.fasta.file=paste0(base.dir,targ, '.fasta')
#    #writeXStringSet(target.fasta, target.fasta.file)
#
#   for(samp in names(targets)) {
#    #samp=targets[[targ]]$samps[1]
#    target.fasta=targets[[samp]]$target.fasta
#    target.fasta.file=paste0(base.dir,samp, '.fasta')
#    writeXStringSet(target.fasta, target.fasta.file)
#
#    pear.file        = paste0(base.dir, samp, '_sort.txt')
#    pear.fasta.file   =paste0(pear.file, '.fasta')
#    blat.table.file   =paste0(pear.file, '.blat.psl')
#    # convert pear output to fasta file
#    convertPear2FASTA(pear.file, pear.fasta.file)
#
#    # run blat with default psl output 
#    runBlat(target.fasta.file=target.fasta.file, pear.fasta.file=pear.fasta.file, blat.table.file=blat.table.file)
#
#    # read in sequences
#    seq=readDNAStringSet(pear.fasta.file)
#    # read in blat alignment stats
#    blat.table=readBlat(blat.table.file)
#   
#    # find reads without indels
#    noIndels=blat.table$tNumInsert==0 & blat.table$qNumInsert==0
#    Indels=!noIndels
#    blat.snps=blat.table[noIndels,]
#
#    #print(targ)
#    print(samp)
#    print('total seqs:')
#    total.seqs=sum(blat.table$seq.counts)
#    print(total.seqs)
#    total.indels=sum(blat.table$seq.counts[Indels])
#    print('total match/mismatch reads:')
#    print(total.seqs-total.indels)
#    print('total indel reads:')
#    print(total.indels)
#    print('fraction of indel reads:')
#    print(total.indels/total.seqs)
#
#    # clip sequences to bits that aligned
#    clipped_seq=padAndClip(seq[as.character(blat.snps$qName)], 
#                   IRanges(start=blat.snps$qStart+1, end=blat.snps$qEnd+1),
#                   Lpadding.letter="N", Rpadding.letter="N")
#    # and align them
#    stacked_seq=stackStrings(clipped_seq, 1,nchar(target.fasta), shift=blat.snps$tStart,
#                             Lpadding.letter="N", Rpadding.letter="N")
#    # restore to original counts
#    stacked_seq=rep(stacked_seq,blat.snps$seq.counts)
#
#    # count up bases at each position
#    snpConsensus=consensusMatrix(stacked_seq, baseOnly=TRUE)
#    rownames(snpConsensus)[5]='N'
#
#    ref.sequence=s2c(as.character(target.fasta[[1]]))
#
#    sname=paste0(targets[[samp]]$info, collapse='_')
#    consensus.match=snpConsensus[cbind(match(ref.sequence, rownames(snpConsensus)), 1:length(ref.sequence))]/apply(snpConsensus,2,sum)
#    plot(consensus.match, ylab='fraction of bases matching target reference', ylim=c(0,1.05), main=sname)
#    text(1:length(ref.sequence), 1.05, ref.sequence, cex=.5)
#
#    #make seqLogo around variant site
#    pwm2=makeSeqLogo(snpConsensus)
#    
#    tinds=strsplit(as.character(blat.table[Indels,]$tStarts),',')
#    tinds=rep(tinds, blat.table$seq.counts[Indels])
#    utinds=as.numeric(unlist(tinds))
#    utinds=utinds[-which(utinds==0)]
#    hist(utinds, breaks=1000, main=paste(sname, 'indel locations') )
#    targets[[samp]]$pear.file        = pear.file        
#    targets[[samp]]$pear.fasta.file  = pear.fasta.file   
#    targets[[samp]]$blat.table.file  = blat.table.file  
#    targets[[samp]]$blat.table       = blat.table       
#    targets[[samp]]$snpConsensus     = snpConsensus     
#    targets[[samp]]$pwm2             = pwm2      
#    targets[[samp]]$utinds           =sort(utinds)    
#
#   }
#dev.off()
#
#
##fraction matching
#results=list()
#for(samp in names(targets)) {
#    ref.sequence=s2c(as.character(targets[[samp]]$target.fasta))
#    snpConsensus=targets[[samp]]$snpConsensus 
#    print(samp)
#    consensus.inds=cbind(match(ref.sequence, rownames(snpConsensus)), 1:length(ref.sequence))
#    consensus.match=snpConsensus[cbind(match(ref.sequence, rownames(snpConsensus)), 1:length(ref.sequence))]/apply(snpConsensus,2,sum)
#    print(targets[[samp]]$info$start+which(consensus.match<.1))
#    print(1-consensus.match[which(consensus.match<.1)])
#    print(snpConsensus[,consensus.inds[,2]][,which(consensus.match<.1)])
#    flag=(consensus.match<.1)+0
#    cons.colname=paste(targets[[samp]]$info$chr, targets[[samp]]$info$start:targets[[samp]]$info$end, ref.sequence, flag, sep=':')
#    consMatrix=snpConsensus
#    colnames(consMatrix)=cons.colname
#    consMatrix2=rbind(consMatrix, rep(0, ncol(consMatrix)))
#    rownames(consMatrix2)[6]='Indels'
#    utinds=targets[[samp]]$utinds
#    tabin=table(cons.colname[utinds])
#    consMatrix2[6, match(names(tabin), colnames(consMatrix))]= as.vector(tabin)
#    results[[samp]]=consMatrix2
#}
#
#
#
#
#save(results, file='/data/coupled_CRISPR/CoupledTest2/consensusMatrices.RData')
#
##pdf(file='~/Desktop/coupled_crispr_edit_distance_test.pdf', width=11, height=8)
##
##
##colfunc <- colorRampPalette(c("red", "purple"), bias=.2)
##colsL=colfunc(4)
##targ='B'
##target.fasta=targets$target.fasta
##ref.sequence=s2c(as.character(target.fasta[[1]]))
##snpConsensus=targets[[3]]$snpConsensus 
##
##samp=(targets$samps)
##consensus.match=snpConsensus[cbind(match(ref.sequence, rownames(snpConsensus)), 1:length(ref.sequence))]/apply(snpConsensus,2,sum)
##max.seq=rownames(snpConsensus)[apply(snpConsensus, 2, which.max)]
##par(yaxs='i')
##plot((1:length(consensus.match))-135, 1-consensus.match, 
##     ylab='fraction of bases', ylim=c(-.05,1.00), main=names(targets$target.fasta),
##         xlim=(c(-18,18)),pch=max.seq, cex=3, lwd=2, xlab='dist from PAM edit'
##         )
##ii=1
##for(n in 4:6){
##    target.fasta=targets$target.fasta
##    ref.sequence=s2c(as.character(target.fasta[[1]]))
##    snpConsensus=targets[[n]]$snpConsensus 
##    max.seq=rownames(snpConsensus)[apply(snpConsensus, 2, which.max)]
##    samp=(targets$samps)[n]
##    consensus.match=snpConsensus[cbind(match(ref.sequence, rownames(snpConsensus)), 1:length(ref.sequence))]/apply(snpConsensus,2,sum)
##    points((1:length(consensus.match))-135, 1-consensus.match, ylab='fraction of bases matching target reference', ylim=c(0,1.05), main=paste(targ),
##         xlim=(c(-20,20)), pch=max.seq, col=colsL[ii], cex=3
##         )
##    ii=ii+1
##    print(n)
##    #readline()
##}
##
##
##
##
##
##ii=1
##for(n in 5:9){
##    target.fasta=targets[[targ]]$target.fasta
##    ref.sequence=s2c(as.character(target.fasta[[1]]))
##    snpConsensus=targets[[targ]][[n]]$snpConsensus 
##    max.seq=rownames(snpConsensus)[apply(snpConsensus, 2, which.max)]
##    samp=names(targets[[targ]])[n]
##    consensus.match=snpConsensus[cbind(match(ref.sequence, rownames(snpConsensus)), 1:length(ref.sequence))]/apply(snpConsensus,2,sum)
##    points((1:length(consensus.match))-135, 1-consensus.match, ylab='fraction of bases matching target reference', ylim=c(0,1.1), main=paste(targ),
##         xlim=(c(-20,20)), pch=max.seq, col=colsR[ii], cex=3
##         )
##    ii=ii+1
##    #print(n)
##    #readline()
##}
##
##dev.off()
##
##
##
##
##
##
##
##for(n in 4:length(names(targets[[targ]]))){
##    target.fasta=targets[[targ]]$target.fasta
##    ref.sequence=s2c(as.character(target.fasta[[1]]))
##    snpConsensus=targets[[targ]][[n]]$snpConsensus 
##
##    samp=names(targets[[targ]])[n]
##    consensus.match=snpConsensus[cbind(match(ref.sequence, rownames(snpConsensus)), 1:length(ref.sequence))]/apply(snpConsensus,2,sum)
##    points((1:length(consensus.match))-135, 1-consensus.match, ylab='fraction of bases matching target reference', ylim=c(0,1.05), main=paste(targ),
##         xlim=rev(c(-20,20)), pch=n
##         )
##    print(n)
##    readline()
##}
##
##
##
##
##
##
##
###    text((1:length(consensus.match))-149, 1:length(ref.sequence), 1.05, ref.sequence, cex=.5)
## #   n=1
##
##
##
##    # some summary info about indels
##    total.indels=sum(blat.table$seq.counts[Indels])
##    total.seqs=sum(blat.table$seq.counts)
##    print(total.indels/total.seqs)
##
##    #ru=rle(sort(utinds))
###barplot(ru$lengths, names=ru$values)
##
##
##targ='A'
##n=4
###Pam edit only in a different color
##colfunc <- colorRampPalette(c("red", "orange"), bias=.2)
##colsR=colfunc(6)
##colfunc <- colorRampPalette(c("red", "purple"), bias=.2)
##colsL=colfunc(4)
##
##target.fasta=targets[[targ]]$target.fasta
##ref.sequence=s2c(as.character(target.fasta[[1]]))
##snpConsensus=targets[[targ]][[n]]$snpConsensus 
##
##samp=names(targets[[targ]])[n]
##consensus.match=snpConsensus[cbind(match(ref.sequence, rownames(snpConsensus)), 1:length(ref.sequence))]/apply(snpConsensus,2,sum)
##max.seq=rownames(snpConsensus)[apply(snpConsensus, 2, which.max)]
##par(yaxs='i')
##plot((1:length(consensus.match))-149, 1-consensus.match, 
##     ylab='fraction of bases', ylim=c(-.05,1.00), main=names(targets[['A']]$target.fasta),
##         xlim=c(-18,18),pch=max.seq, cex=3, lwd=2, xlab='dist from PAM edit'
##         )
##for(n in 5:10){
##    target.fasta=targets[[targ]]$target.fasta
##    ref.sequence=s2c(as.character(target.fasta[[1]]))
##    snpConsensus=targets[[targ]][[n]]$snpConsensus 
##    max.seq=rownames(snpConsensus)[apply(snpConsensus, 2, which.max)]
##    samp=names(targets[[targ]])[n]
##    consensus.match=snpConsensus[cbind(match(ref.sequence, rownames(snpConsensus)), 1:length(ref.sequence))]/apply(snpConsensus,2,sum)
##    points((1:length(consensus.match))-149, 1-consensus.match, ylab='fraction of bases matching target reference', ylim=c(0,1.05), main=paste(targ),
##         xlim=c(-20,20), pch=max.seq, col=colsR[n-4], cex=3
##         )
##    #print(n)
##   # readline()
##}
##ii=1
##for(n in c(13,14,11,12)){
##    target.fasta=targets[[targ]]$target.fasta
##    ref.sequence=s2c(as.character(target.fasta[[1]]))
##    snpConsensus=targets[[targ]][[n]]$snpConsensus 
##    max.seq=rownames(snpConsensus)[apply(snpConsensus, 2, which.max)]
##    samp=names(targets[[targ]])[n]
##    consensus.match=snpConsensus[cbind(match(ref.sequence, rownames(snpConsensus)), 1:length(ref.sequence))]/apply(snpConsensus,2,sum)
##    points((1:length(consensus.match))-149, 1-consensus.match, ylab='fraction of bases matching target reference', ylim=c(0,1.05), main=paste(targ),
##         xlim=c(-20,20), pch=max.seq, col=colsL[ii], cex=3
##         )
##    ii=ii+1
##    #print(n)
##    #readline()
##}
##
##
##
###Pam edit only in a different color
##colfunc <- colorRampPalette(c("red", "orange"), bias=.2)
##colsR=colfunc(5)
##colfunc <- colorRampPalette(c("red", "purple"), bias=.2)
##colsL=colfunc(4)
#


