# RNA-seq pipeline
#Explanation of every parameter used using the respective tools manual and other external sources

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

#step3 multiQC

#In this pipeline we are mainly using multiqc to interpret the data from cutadapt and fastQC.
#The report produced is in html format with various reports from different bioinformatics tools such as cutadpat and fastQC. Furthermore multiqc produces data in the form of multiqc_data that can be used for downstream analysis and include excel files that contain more information than the html report.
#In order to configure how multiqc works, we can use the command multiqc --help; and to ignore and replace a current report, we can use the flag command multiqc --force . (current directory).
#The reports produced are usually tab separated but we can change the data format to jason using the command: multiqc --data-format  json.
#download matplotlib and test 04 fastqc report that has been downloaded
multiqc fastqc_report(note: folder name)

#step4 sequence alignment (STAR)

#Two main steps for this STAR workfflow which include 1) creating a genome index and 2) mapping reads to th genome.
# --genomeDir command is used to create genome index and provides the path/to/genomeDir (/home-1/bpark27@jhu.edu/sbiswal1/References/mm10). The reference genome is mm10
# --readFilesIn provides the RNA-seq reads (sequences) in the form of fastq files.
# --readFilesCommand takes the file name as input parameter and sends the uncompressed output to stdout. zcat is used here beause files here are gzipped.
# --runThreadN defines the number of threads to be used for genome generation according to the number of available cores (processing unit receiving instructions and performing calculation) on the server node (basic unit)
#--genomeLoad has the default NoShared Memory so that we can avoid using shared memory in super computer.
# --outFileNamePrefix provides all the output files from STAR in one directory
# --outputFilterMultimapNmax 20 is an ENCODE option that provides the max number of multiple alignments allowed for a read. If exceeded the read is considered unmapped.
# --alignSJoverhangMin 8 is another ENCODE option  that provides the minimun overhang (unpaired nucleotides in the end of DNA molecule) for unannotated junctions
# --alignSJDBoverhangMin 1 is also an ENCODE option that provides the minimun overhang for annotated junctions
# --outFilterMismatchNmax 999 is another ENCODE option that provides the maximum number of mismatches per pair, large number switches off this filter.
# --outFilterMismatchNoverReadLmax 0.04 is another ENCODE option that provides the maximun number of mismatches per pair relative to read length: for 2*100b, max number of mismatches is 0.06*200=8 for the paired read
# --alignIntronMin 20 is another ENCODE option that provides the minimum intron length
# --alignIntronMax 1000000 is another ENCODE option that provides the maximun intron length
# --alignMateGaGapMax 1000000 is the final ENCODE option in STAR manual and provides the maximum genomic distance between mates
# --outSAMheaderCommentFile provides a string with the path to the file with @CO (comments) lines of the SAM header; in this case COfile.txt
# outSAMheaderHD provides a string with @HD (header) line of the SAM header; in this case @HD VN:1.4 SO:coordinate
# --outSAMunmapped within provides the output of the unmapped reads within the main SAM file (i.e Aligned.out.sam)
# --outFilterType BySJout reduces the number of "spurious" (fake) junctions
# --outSAMattributes accepts a list of 2-character SMA attributes. NH HI NM MD have standard meaning as defined in the SAM format specifications and AS id the local alignment score (paired for paired end reads).
# --outSAMtype BAM SortedByCoordinate provides the output sorted by coordinate Aligned.sortedByCoord.out.bam file, similar ti samtools sort command.
# --quantMode TranscriptomeSAM outputs alignments translated into transcript coordinates in the Aligned.toTranscriptome.out.bam file (in addition to alignment in genomic coordinates in Aligned.*.sam/bam files). These transcriptomic alignments can be used with various transcript quantification software that require reads to be mapped to transcriptome, such as RSEM or eXpress.
# --sjdbScore is an int (1 in this case) that provides an extra alignment score for alignments that cross database junctions.
# --limitBAMsortRAM is an int > 0 (60000000000 in this case) that provides maximum available RAM for sorting BAM. if =0, it will be set to the genome index size

STAR --genomeDir /home-1/bpark27@jhu.edu/sbiswal1/References/mm10/ --readFilesIn read1_trimmed.fq read2_trimmed.fq --readFilesCommand zcat --runThreadN 2 --genomeLoad NoSharedMemory -outFileNamePrefix prefix_ --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMheaderCommentFile COfile.txt --outSAMheaderHD @HD VN:1.4 SO:coordinate --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --sjdbScore 1 --limitBAMsortRAM 60000000000

#step5 featureCount

#featureCounts is a general-purpose read summarization function, which assigns mapped reads (RNA-seq reads in this case) to genomic features or meta-features.
#A feature is an interval (range of positions) on one of the reference sequences eg exons
#A meta-feature is a set of features that represents a biological construct of interest eg genes.
#featureCounts can summarize reads at either the feature or meta-feature levels.
# The data input should consist of one or more aligned reads (short or long reads) in either SAM or BAM format and a list of genomic features (GTF in this case). The format of input reads reads is automatically detected (SAM or BAM).
# -T 2 implies that there are 2 threads and hence there 2 processes to be executed concurrently by featureCount
# -s int (0 in this case) implies that in our case read counting is not strand specific (unstranded).' -s 1' implies stranded and'-s 2' implies reversely stranded)
# -t string (in this case exon which is default) is an inbuilt option in featureCount that specifies the feature type. Only rows which have the matched feature type in the provided GFT annotation file will be included for read counting.
# -g string (gene_id) is another default function in featureCount that specifies the attribute used to to group features in the GTF file
# -a string (annotation file) is also an inbuilt option in featureCount wherein in this case the annotation is in gtf format and is from gencode.vM13. Note that when using featureCount the annotation file should be either in GTF, GFF or SAF format.
# -o string is the paramter that indicates that whatever comes after it is the output file that contains the number of reads assigned to each meta-feature.
# We are summarizing BAM format (prefix_Aligned.sortedByCoord.out.bam) unstranded reads (-s 0) using 2 threads (-T 2). 

featureCounts -T 2 -s 0 -t exon -g gene_id -a gencode.vM13.annotation.gtf -o output_count.txt prefix_Aligned.sortedByCoord.out.bam
