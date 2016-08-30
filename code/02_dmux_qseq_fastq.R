### Core processing code to convert Hiseq fastq into alignment files ===============================================
#
# 1) demultiplex fastq.gz files obtained from BSCRC Hiseq using 1 or 2 index reads
# 2) run PEAR to assemble overlapping reads
# 3) trim 'adapters' and constant sequences in amplicon library
# 4) convert to PHRED33 (perhaps this would be easier during the initial conversion from qseq to fastq)
# 5) extract internal barcodes, classify reads by internal barcode, and add this information to the read name
# 6) align to expected repair template sequence
#===================================================================================================================

library(ShortRead)
library(seqinr)
library(doMC)
n.cores=33
#n.cores=1
registerDoMC(cores=n.cores)
cl <- makeCluster(getOption("cl.cores", n.cores))

#batch='eQTL3' ;
lane.identifier='082316_Estops_Hiseq' #042916_Estops_Hiseq'

readformat= 'gz.fastq'
base.dir='/media/jbloom/d1/coupled_CRISPR/Experiments/082316/'
out.base.dir = base.dir

sample.sheet.file = paste0(base.dir, 'Reference/2016.08.06_MS_TruSeq_hiseq.csv')
#q/qmedia/jbloom/d1/coupled_CRISPR/042916_Estops_Hiseq/2016.04.15_MS_tag_v_pcr.csv'

base.dir.qseq=base.dir #(paste(base.dir, 'qseq/', batch, '/', sep=''))

#read in a miseq formatted sample sheet
ss=read.delim(sample.sheet.file, header=T, sep=',', stringsAsFactors=F, skip=19) #skip=20 # if miseq sample sheet

# additional to code to add reverse complement to sample sheet 
ss$index=toupper(ss$index)
# ss2=ss
# ss2$index=as.character(reverseComplement(DNAStringSet(ss2$index)))
# ss2$Sample_Name=paste0(ss2$Sample_Name, '_rc')
# ss=rbind(ss,ss2)
#--------------------------------------------------------------------------------------------

n700.inds=unique(ss$index)
n500.inds=unique(ss$index2)

ind384=paste(match(ss$index, n700.inds)) #, match(ss$index2, n500.inds), sep='_')

#dir.create(paste0(out.base.dir, '/fastqrc/'))
#dir.create(paste0(out.base.dir, 'fastqrc/Indices/'))

dir.create(paste0(out.base.dir, '/fastq/'))
dir.create(paste0(out.base.dir, 'fastq/Indices/'))
#additional key
#key2=read.delim(file='/media/jbloom/d1/coupled_CRISPR/042916_Estops_Hiseq/sample_key.txt', header=F, sep='\t')
#key2.names=apply(key2[,2:4], 1, paste, collapse='_')
#key2.names=gsub('\\s','', key2.names)
#ss$Sample_Plate=key2.names[ss$Sample_Name]


# Version starting with gzipped fastq files  -------------------------------------------------------------------------
i1=paste(base.dir.qseq, lane.identifier, '_index1.fastq.gz', sep='')
#i2=paste(base.dir.qseq, lane.identifier, '_index2.fastq.gz', sep='')
r1=paste(base.dir.qseq, lane.identifier, '_read1.fastq.gz', sep='')
r2=paste(base.dir.qseq, lane.identifier, '_read2.fastq.gz', sep='')


nbuffer=2e6

#rf1=rle(sort(as.character(sread(fq.i1))))
#rf2=rle(sort(as.character(sread(fq.i2))))

#o1=rf1$values[order(rf1$lengths, decreasing=T)]
#o1=data.frame(o1, rf1$lengths[order(rf1$lengths, decreasing=T)])

#o2=rf2$values[order(rf2$lengths, decreasing=T)]
#o2=data.frame(o2, rf2$lengths[order(rf2$lengths, decreasing=T)])

#PRELOADED to here
fi1 =FastqStreamer(i1, nbuffer, readerBlockSize=1e7,verbose=F)
#fi2 =FastqStreamer(i2, nbuffer, readerBlockSize=1e7)
fr1 =FastqStreamer(r1, nbuffer, readerBlockSize=1e7)
fr2 =FastqStreamer(r2, nbuffer, readerBlockSize=1e7)
repeat {
    fq.i1=yield(fi1)   
    #fq.i2=yield(fi2)   
    fq.r1=yield(fr1)   
    fq.r2=yield(fr2)   

    if(length(fq.i1) ==0 ) {break }
    print('calculating distance to index1')
    # calculate distance matrix of reads indices (read1 is N700 read)
    index1.dist.calc=srdistance(sread(fq.i1), n700.inds	)
    d1mat=do.call('cbind', index1.dist.calc)
    # find closest match for each read
    d1min.ind=parRapply(cl, d1mat, which.min)
    # find edit distance of closest match
    d1min=d1mat[cbind(1:nrow(d1mat), d1min.ind)]

#     # lop off the last base
#     # s7=DNAStringSet(sread(fq.i2),start=1, end=8)
#     print('calculating distance to index2')
#     # calculate distance matrix of reads indices (read2 is N500 read)
#     index2.dist.calc=srdistance(sread(fq.i2), n500.inds)
#     d2mat=do.call('cbind', index2.dist.calc) 
#     # find closest match for each read
#     d2min.ind=parRapply(cl, d2mat, which.min)
#     # find edit distance of closest match
#     d2min=d2mat[cbind(1:nrow(d2mat), d2min.ind)]

    match.indices=d1min<2
    keep.reads=which(match.indices) #& d2min<2)
    print(paste( (length(keep.reads)/nbuffer)*100 ,  ' % of reads retained'))
    d1.split.ind=d1min.ind[keep.reads]
    #d2.split.ind=d2min.ind[keep.reads]
    matchtoinput=paste(d1.split.ind ) #, d2.split.ind, sep='_')
    
    #    readID=paste(ss$Sample_Name, ss$Sample_Plate, sep='-')[match(matchtoinput, ind384)]
    readID=paste(ss$Sample_Name)[match(matchtoinput, ind384)]
    readID.split=split(keep.reads, readID)

    unmatched.reads=which(!match.indices)
    fread1.unmatched=paste(out.base.dir, 'fastq/', 'unmatched', '_R1.fq.gz', sep='')
    fread2.unmatched=paste(out.base.dir, 'fastq/', 'unmatched', '_R2.fq.gz', sep='')
    iread1.unmatched=paste(out.base.dir, 'fastq/Indices/', 'unmatched', '_I7.fq.gz', sep='')
    writeFastq(fq.r1[unmatched.reads], file=fread1.unmatched, mode='a', full=FALSE, compress=TRUE)
    writeFastq(fq.r2[unmatched.reads], file=fread2.unmatched, mode='a', full=FALSE, compress=TRUE)
    writeFastq(fq.i1[unmatched.reads], file=iread1.unmatched, mode='a', full=FALSE, compress=TRUE)


    foreach(i=1:length(readID.split)) %dopar% {
        #file1=paste(out.base.dir, 'fastq/', names(readID.split)[i] ,'-', lane.identifier, '.fq.gz', sep='')
        file1=paste(out.base.dir, 'fastq/', names(readID.split)[i], '_R1.fq.gz', sep='')
        file2=paste(out.base.dir, 'fastq/', names(readID.split)[i], '_R2.fq.gz', sep='')
        file3=paste(out.base.dir, 'fastq/Indices/', names(readID.split)[i], ':I7.fq.gz', sep='')
        #file4=paste(out.base.dir, 'fastq/Indices/', names(readID.split)[i], ':I5.fq.gz', sep='')

        writeFastq(fq.r1[readID.split[[i]]], file=file1, mode='a', full=FALSE, compress=TRUE)
        writeFastq(fq.r2[readID.split[[i]]], file=file2, mode='a', full=FALSE, compress=TRUE)
        writeFastq(fq.i1[readID.split[[i]]], file=file3, mode='a', full=FALSE, compress=TRUE)
        #writeFastq(fq.i2[readID.split[[i]]], file=file4, mode='a', full=FALSE, compress=TRUE)
    }
}
#-----------------------------------------------------------------------------------------------------------------------




##### # for miseq Start HERE !!!! ###############################
## No need to demultiplex with miseq ####
#find /home/jbloom/Downloads/miseq/  -type f -iname '*.fastq.gz' -exec cp {} /media/jbloom/d1/coupled_CRISPR/080316_fixed_structural_element/fastq/ \;
#find /home/jbloom/MS_TruSeq-31815832/  -type f -iname '*.fastq.gz' -exec cp {} /media/jbloom/d1/coupled_CRISPR/080816_truseq/fastq/ \;

#There's only one difference from last time: samples 13 and 14 switch places with 15 and 16, respectively.

# run pear to resolve overlapping paired reads into one read ------------------------------------------------------------
# note this assumes illumina 1.3 Quality encoding (64bit) (BSCRC Hiseq 4000)

#readID=paste(ss$Sample_Name, ss$Sample_Plate, sep='-')
out.base.dir='/media/jbloom/d1/coupled_CRISPR/Experiments/082316/'
readID=list.files('/media/jbloom/d1/coupled_CRISPR/Experiments/082316/fastq/')
readID=unique(gsub('_R.*.fq.gz', '', readID))
readID=readID[-grep('Indices', readID)]
#c(seq(1,16) ) #, 'plate_gal','plate_glu')
dir.create(paste0(out.base.dir, 'pear/'))
for (rr in readID) {
    ffq=paste(out.base.dir, 'fastq/', rr, '_R1.fq.gz', sep='')
    rfq=paste(out.base.dir, 'fastq/', rr, '_R2.fq.gz', sep='')
    ofq=paste(out.base.dir, 'pear/', rr, '_pear_out', sep='')
    #system(paste('pear -f', ffq , '-r' , rfq, '-o' , ofq , '-v 20 -m 2750 -j 10' )) # -b 64' ) )
    # added --phred-base 64 for Illumina 1.3 format 
    system(paste('/home/jbloom/Local/bin/pear -f', ffq , '-r' , rfq, '-o' , ofq , '-v 20 -m 275 -j 70 -b 64' ))
}
#-----------------------------------------------------------------------------------------------------------------------

# remove common adapter sequences --------------------------------------------------------------------------------------

seqL1='CCTTCTTGACGAGTTCTTCTGAATTATTAACGCTTACAATTTCCTGATGCGGTATTTTCTCCTTACGCATCTGTGCGGTATTTCACACCGCTCGGAGGGCAACTATCC'
seqR1='CCAGCTTCACACGGC'
#seqL1='GCGGTATTTCACACCGCATACGTCAGATGTGTATAAGAGACAG'
# for 062916
#seqL1='CCTTCTTGACGAGTTCTTCTGAATTATTAACGCTTACAATTTCCTGATGCGGTATTTTCTCCTTACGCATCTGTGCGGTATTTCACACCGCTCGGAGGGCAACTATCC'
# for 080316
# followed by 12 bp barcode
# then 6 bp fixed sequence ACGCGT
# then one of these 
#seqR1='CCAGCTTCACACGGCCGGCCGGTACCCAATTCGCCCTATAGTGAGTCGTAT'
#seqR2='CCAGCTTCACACGGCCGGTACCCAATTCGCCCTATAGTGAGTCGTAT'
# for 080316
#seqL1='TCGGAGGGCAACTATCC'
#seqR1='CCAGCTTCACACGGC'
foreach(rr = readID)  %dopar% {
    fin=paste(out.base.dir, 'pear/', rr, '_pear_out.assembled.fastq', sep='')
    fout=paste(out.base.dir, 'pear/', rr, '_adapter_trimmed.fastq', sep='')
    system(paste('~/.local/bin/cutadapt  -n 3 -g', seqL1, '-b', seqR1, '-m 40 -o', fout, fin))
                 ###'-b', seqR2, '-m 40 -o', fout, fin))

}
#-----------------------------------------------------------------------------------------------------------------------

# convert to phred33 and compress--------------------------------------------------------------------------------------
foreach(rr = readID) %dopar% {
    fin = paste(out.base.dir, 'pear/', rr, '_adapter_trimmed.fastq', sep='')
    fout = paste(out.base.dir, 'pear/', rr, '_adapter_trimmed_phred33.fq.gz', sep='')
    system(paste('/home/jbloom/Local/seqtk/seqtk seq -Q64 -V', fin, '| pigz - > ', fout))
    #system(paste('pigz', fout ) )#, ' > ', fout))
}
#----------------------------------------------------------------------------------------------------------------------


swws=(cbind(amb('S'), amb('W'), amb('W'), amb('S')))
toupper(swws)
eg=expand.grid(c(1,2),c(1,2),c(1,2),c(1,2))
swws=toupper(apply(t(apply(eg, 1, function(x) swws[cbind(c(x),c(1,2,3,4))])), 1, c2s))

wssw=(cbind(amb('W'), amb('S'), amb('S'), amb('W')))
toupper(wssw)
eg=expand.grid(c(1,2),c(1,2),c(1,2),c(1,2))
wssw=toupper(apply(t(apply(eg, 1, function(x) wssw[cbind(c(x),c(1,2,3,4))])), 1, c2s))


# Now trim off internal barcode and attach to readname ----------------------------------------------------------------
# Also, trim off fixed sequence after barcode 
# create new directory structure for output
dir.create(paste0(out.base.dir, 'processed/'))
dir.create(paste0(out.base.dir, 'processed/fastq/'))
dir.create(paste0(out.base.dir, 'processed/bam/'))
#NNNNNNNNSWWS
#NNNNNNNNSWWS
#ACGCGT
barcode.length = 12 #5
fixed.sequence.length=6
nbuffer=5e6
for(rr in readID) {
    #PRELOADED to here
    print(rr)
    #rr.file=paste(out.base.dir, 'pear/', rr, 'CCTGTTTGTATGCAATTTTTAATCAGTTTTCCTCTCTGACACCTTTCTATCACTTTTTTTTTCATATCTGCCCTCAGGTCGCATCCTCTTTCTTCCCTACC_adapter_trimmed_phred33.fq.gz', sep='')
    rr.file=paste(out.base.dir, 'pear/', rr, '_adapter_trimmed_phred33.fq.gz', sep='')

    # this generalizes to very large fastq files 
    fi1 =FastqStreamer(rr.file, nbuffer, readerBlockSize=1e7,verbose=F)
    repeat {
        #ShortReadQ class 
        rfq=yield(fi1)   
        if(length(rfq) ==0 ) {break }
        # read the whole thing into memory ... this can obviously go very badly but should provide a 
        #template for a more memory efficient version     #rfq=readFastq(rr.file, withIDs=TRUE)
        # ---------------------------------------------------------------------------------------------
        cread=sread(rfq)
        barcodes=subseq(cread, start=1, width=barcode.length) #5)
        barcode.expt=subseq(barcodes, start=9, width=4)
        barcode.id=rep('-', length(barcode.expt))
        barcode.swws = as.character(barcode.expt) %in% swws
        barcode.wssw = as.character(barcode.expt) %in% wssw
        barcode.id[barcode.swws]='swws'
        barcode.id[barcode.wssw]='wssw'


        trimmed.read=subseq(cread, start=barcode.length+fixed.sequence.length+1) #19) #start=12 )
        # fixed sequence ACGCGT
        cname=as.character(id(rfq))
        cname=gsub('\\s', ':', cname)
        cname=paste(cname, barcodes, barcode.id, sep=':')
        newname=BStringSet(cname)
        rfout=paste(out.base.dir, 'processed/fastq/', rr, '_adapter_trimmed_phred33_barcode.fq.gz', sep='')
        writeFastq( ShortReadQ(sread=trimmed.read, quality=narrow(quality(rfq), start=barcode.length+fixed.sequence.length+1), id= newname), 
              file=rfout, mode='a', full=FALSE, compress=TRUE)
    }
    close(fi1)
}
#----------------------------------------------------------------------------------------------------------------------


# mapping -----------------------------------------------------------------------------------------------------------------
#file was here /media/jbloom/d1/coupled_CRISPR/062916_Reworked/ref/eStops_repair_template_only.fasta'
reference.sequences='/media/jbloom/d1/coupled_CRISPR/Reference/eStops_repair_template_only.fasta'
for( rr in readID) {
        rfin=paste(out.base.dir, 'processed/fastq/', rr, '_adapter_trimmed_phred33_barcode.fq.gz', sep='')
        rout=paste(out.base.dir, 'processed/bam/', rr, sep='')
        system(paste('bwa mem -t 70', reference.sequences, rfin , ' | samtools view -hSu - | samtools sort - ', rout))
}
#----------------------------------------------------------------------------------------------------------------------------


# index bam files  -----------------------------------------------------------------------------------------------------
for( rr in readID) {
    bam.file=paste(out.base.dir, 'processed/bam/', rr, '.bam', sep='')
    system(paste('samtools index', bam.file))
 }
# --------------------------------------------------------------------------------------------------------------------
