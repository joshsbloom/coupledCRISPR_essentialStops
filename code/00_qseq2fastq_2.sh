## Pull data from pan.pelligrini.mcdb.ucla.edu

# first sequencing
##SxaQSEQsXap117L1:QF3Gf9kB8xf3 
##SxaQSEQsXap117L2:mn8fq6fK8PA3

# resequence
#SxaQSEQsXbp082L1:Qr8rj4CV4uW8
#SxaQSEQsXbp082L2:jv7YQ4HN3Fq4

#rsync --recursive --times --verbose --stats --progress --itemize-changes --sockopts=SO_RCVBUF=3145728,SO_SNDBUF=3145728 rsync://SxaQSEQsXap117L1@pan.pellegrini.mcdb.ucla.edu/SxaQSEQsXap117L1 L1/ 
#rsync --recursive --times --verbose --stats --progress --itemize-changes --sockopts=SO_RCVBUF=3145728,SO_SNDBUF=3145728 rsync://SxaQSEQsXap117L2@pan.pellegrini.mcdb.ucla.edu/SxaQSEQsXap117L2 L2/ 

#rsync --recursive --times --verbose --stats --progress --itemize-changes --sockopts=SO_RCVBUF=3145728,SO_SNDBUF=3145728 rsync://SxaQSEQsXbp082L1@pan.pellegrini.mcdb.ucla.edu/SxaQSEQsXbp082L1 L1/ 
#rsync --recursive --times --verbose --stats --progress --itemize-changes --sockopts=SO_RCVBUF=3145728,SO_SNDBUF=3145728 rsync://SxaQSEQsXbp082L2@pan.pellegrini.mcdb.ucla.edu/SxaQSEQsXbp082L2 L2/ 


#### to extract just the sequence data
###awk 'NR%4==2' a

## replicate hiseq runs, combine the two runs together and use multicore gzip implementation to compress
#cat L1_index1.fastq L2_index1.fastq | pigz - > 082316_Estops_Hiseq_index1.fastq.gz
#cat L1_read1.fastq L2_read1.fastq   | pigz - > 082316_Estops_Hiseq_read1.fastq.gz
#cat L1_read2.fastq L2_read2.fastq   | pigz - > 082316_Estops_Hiseq_read2.fastq.gz

cat L1_2_index1.fastq L2_2_index1.fastq | pigz - > 090916_Estops_Hiseq_index1.fastq.gz
cat L1_2_read1.fastq L2_2_read1.fastq   | pigz - > 090916_Estops_Hiseq_read1.fastq.gz
cat L1_2_read2.fastq L2_2_read2.fastq   | pigz - > 090916_Estops_Hiseq_read2.fastq.gz

## see code/SJC-to-UCLA-HiSeq-Users-20120125-V2.pdf for the logic of the BASH script for converting qseq to fastq
INdir='/media/jbloom/d1/coupled_CRISPR/Experiments/082316/L1/'
IN1='1_1'
IN2='1_2'
IN3='1_3'
#IN4='1_4'

OUTfile1='/media/jbloom/d1/coupled_CRISPR/Experiments/082316/L1_2_read1.fastq'
OUTfile2='/media/jbloom/d1/coupled_CRISPR/Experiments/082316/L1_2_index1.fastq'
OUTfile3='/media/jbloom/d1/coupled_CRISPR/Experiments/082316/L1_2_read2.fastq'
#OUTfile4='/media/jbloom/d1/coupled_CRISPR/042916_Estops_Hiseq/L1_read2.fastq'

# note filter for reads_passed_filter is $11==1

(for IN in "${INdir}"/s_"${IN1}"*.txt.gz; do
   (echo '@@@' "Working on ${IN}...  $(date)") > /dev/stderr
   gunzip -c  "${IN}" | awk -F $'\t' 'BEGIN { serial=0; }           \
     { serial++; gsub(/\./,"N", $9);                                \
      if ($11==1) printf("@%s_%s_%s_%s_%09d_%s\n%s\n+\n%s\n",       \
       $2,$3,$4,$11,serial,$8, $9, $10); }'; done) > "${OUTfile1}" 
(echo '@@@' "Done.  $(date)") > /dev/stderr    

(for IN in "${INdir}"/s_"${IN2}"*.txt.gz; do
   (echo '@@@' "Working on ${IN}...  $(date)") > /dev/stderr
   gunzip -c  "${IN}" | awk -F $'\t' 'BEGIN { serial=0; }           \
     { serial++; gsub(/\./,"N", $9);                                \
      if ($11==1) printf("@%s_%s_%s_%s_%09d_%s\n%s\n+\n%s\n",       \
       $2,$3,$4,$11,serial,$8, $9, $10); }'; done) > "${OUTfile2}" 
(echo '@@@' "Done.  $(date)") > /dev/stderr  

(for IN in "${INdir}"/s_"${IN3}"*.txt.gz; do
   (echo '@@@' "Working on ${IN}...  $(date)") > /dev/stderr
   gunzip -c  "${IN}" | awk -F $'\t' 'BEGIN { serial=0; }           \
     { serial++; gsub(/\./,"N", $9);                                \
      if ($11==1) printf("@%s_%s_%s_%s_%09d_%s\n%s\n+\n%s\n",       \
       $2,$3,$4,$11,serial,$8, $9, $10); }'; done) > "${OUTfile3}" 
(echo '@@@' "Done.  $(date)") > /dev/stderr  

(for IN in "${INdir}"/s_"${IN4}"*.txt.gz; do
   (echo '@@@' "Working on ${IN}...  $(date)") > /dev/stderr
   gunzip -c  "${IN}" | awk -F $'\t' 'BEGIN { serial=0; }           \
     { serial++; gsub(/\./,"N", $9);                                \
      if ($11==1) printf("@%s_%s_%s_%s_%09d_%s\n%s\n+\n%s\n",       \
       $2,$3,$4,$11,serial,$8, $9, $10); }'; done) > "${OUTfile4}" 
(echo '@@@' "Done.  $(date)") > /dev/stderr  



INdir='/media/jbloom/d1/coupled_CRISPR/Experiments/082316/L2/'
IN1='2_1'
IN2='2_2'
IN3='2_3'
#IN4='1_4'

OUTfile1='/media/jbloom/d1/coupled_CRISPR/Experiments/082316/L2_2_read1.fastq'
OUTfile2='/media/jbloom/d1/coupled_CRISPR/Experiments/082316/L2_2_index1.fastq'
OUTfile3='/media/jbloom/d1/coupled_CRISPR/Experiments/082316/L2_2_read2.fastq'
#OUTfile4='/media/jbloom/d1/coupled_CRISPR/042916_Estops_Hiseq/L1_read2.fastq'

# note filter for reads_passed_filter is $11==1

(for IN in "${INdir}"/s_"${IN1}"*.txt.gz; do
   (echo '@@@' "Working on ${IN}...  $(date)") > /dev/stderr
   gunzip -c  "${IN}" | awk -F $'\t' 'BEGIN { serial=0; }           \
     { serial++; gsub(/\./,"N", $9);                                \
      if ($11==1) printf("@%s_%s_%s_%s_%09d_%s\n%s\n+\n%s\n",       \
       $2,$3,$4,$11,serial,$8, $9, $10); }'; done) > "${OUTfile1}" 
(echo '@@@' "Done.  $(date)") > /dev/stderr    

(for IN in "${INdir}"/s_"${IN2}"*.txt.gz; do
   (echo '@@@' "Working on ${IN}...  $(date)") > /dev/stderr
   gunzip -c  "${IN}" | awk -F $'\t' 'BEGIN { serial=0; }           \
     { serial++; gsub(/\./,"N", $9);                                \
      if ($11==1) printf("@%s_%s_%s_%s_%09d_%s\n%s\n+\n%s\n",       \
       $2,$3,$4,$11,serial,$8, $9, $10); }'; done) > "${OUTfile2}" 
(echo '@@@' "Done.  $(date)") > /dev/stderr  

(for IN in "${INdir}"/s_"${IN3}"*.txt.gz; do
   (echo '@@@' "Working on ${IN}...  $(date)") > /dev/stderr
   gunzip -c  "${IN}" | awk -F $'\t' 'BEGIN { serial=0; }           \
     { serial++; gsub(/\./,"N", $9);                                \
      if ($11==1) printf("@%s_%s_%s_%s_%09d_%s\n%s\n+\n%s\n",       \
       $2,$3,$4,$11,serial,$8, $9, $10); }'; done) > "${OUTfile3}" 
(echo '@@@' "Done.  $(date)") > /dev/stderr  

(for IN in "${INdir}"/s_"${IN4}"*.txt.gz; do
   (echo '@@@' "Working on ${IN}...  $(date)") > /dev/stderr
   gunzip -c  "${IN}" | awk -F $'\t' 'BEGIN { serial=0; }           \
     { serial++; gsub(/\./,"N", $9);                                \
      if ($11==1) printf("@%s_%s_%s_%s_%09d_%s\n%s\n+\n%s\n",       \
       $2,$3,$4,$11,serial,$8, $9, $10); }'; done) > "${OUTfile4}" 
(echo '@@@' "Done.  $(date)") > /dev/stderr  

