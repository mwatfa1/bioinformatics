# RNA-seq pipeline

# Step #1
# Trimming raw fastq files using Trim Galore
# What is trim galore?
#Trim galore is a wrapper tool around Cutadapt and FastQC to consistently apply quality and adapter trimming to FastQ files, with some extra functionality for MspI-digested RRBS-type (Reduced Representation Bisufite-Seq) libraries.
# A wrapper's job is to pass along its inputs to some other tool, with slight modifications
#Reduced representation bisulfite sequencing (RRBS) is an efficient and high-throughput technique used to analyze the genome-wide methylation profiles on a single nucleotide level. This technique combines restriction enzymes and bisulfite sequencing in order to enrich for the areas of the genome that have a high CpG content. This technique is used in order to reduce the amount of nucleotides needed to be sequenced to 1% of the genome.

#The parameters used in trim galore are -q and --paired
# -q : Trims low-quality ends from reads. In this case anything lower than 25 is trimmed
# --paired : This option performs length trimming of quality/adapter/RRBS trimmed reads for paired-end files

#fastq files used
#read1.fq = /home-3/mwatfa1@jhu.edu/sbiswal1/RNAseq_Musa/03_trimmed_fastq_files/read1.fq
#read2.fq = /home-3/mwatfa1@jhu.edu/sbiswal1/RNAseq_Musa/03_trimmed_fastq_files/read2.fq
trim_galore -q 25 --paired read1.fq  read2.fq

#step2 fastqc report of trimmed fastq headata (using trim_galore)
#The parameter used here is -o which means that the output of the fastqc report for read1.fq/read2.fq should be in the parent directory of 03_trimmed_fastqc_files/read1.fq and 03_trimmed_fastqc_files/read2.fq respectively.
fastqc -o read2.fq.gz_fastqc ../03_trimmed_fastq_files/read1.fq
fastqc -o read2.fq.gz_fastqc ../03_trimmed_fastq_files/read2.fq
