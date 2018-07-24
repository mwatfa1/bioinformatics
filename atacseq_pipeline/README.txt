This is ATAC-seq pipeline developed by Bongsoo and Musa
2018/07/16

#First check module of bwa available on server;
#1 command: module avail bwa
#module found is : bwa/0.7.17

#Next load bwa version present:
#2 command: module load bwa/0.7.17

#Now run bwa
#3 command: bwa

#To run bwa on server, we need to enable the marcc interact mode
#4 command: salloc -J interact -N 1-1 -n 1 -c 4 --mem=20G --time=360 -p shared srun --pty bash

#make a directory for the reference mm10 and go into the directory
#5 command: mkdir mm10
#6 command: cd mm10

#copy reference genome from server into the new mm10 directory
#7 command: cp ~/sbiswal1/References/mm10/mm10.fa* .

#we can now move up to the previous directory in order to index the r1 and r2 fastq files to the reference genome and output the result onto mapping.sam
#8 command cd ..
#9 command bwa mem -t 4 mm10/mm10.fa r1.fastq r2.fastq > mapping.sam

#load samtools in order to be able to carry out further analysis
#10 command: module load samtools

#convert mapping.sam output to mapping.bam in order to be able to continue the analysis 
#11 command: samtools view -S -b mapping.sam > mapping.bam

#find uniquely mapped reads using samtools
#12 command: samtools view -bq 1 mapping.bam > unique.bam

#remove duplicated reads and create 
#13 command: rmdup unique.bam dup.bam

#The purpose of using MAPQ score is to quantify the probability of correctly mapping a random read, thereby quantifying mapping quality.
#mapq score of 20 = 0.99
#mapq score of 30 = 0.999

