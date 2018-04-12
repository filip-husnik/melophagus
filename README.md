# Commands to reproduce the RNA-Seq analysis of M. ovinus published in Husnik et al.

# Melophagus ovinus
*Data quality assesment in RseQC and FastQC*

*Adapter and quality trimming in Cutadapt and Sickle*

*SSU contamination assesment in PhyloFlash*

```
#RNA-Seq
/opt/phyloFlash_v3.1b2/phyloFlash.pl -lib female -CPUs 12 -read1 ../1-MO3_130830_L004_R0.fastq.gz -read2 ../1-MO3_130830_L004_R1.fastq.gz -dbhome /Data/filip/phyloFlash_DB/128/ -zip -readlength 100 -almosteverything -id 65 -readlimit 5000000
/opt/phyloFlash_v3.1b2/phyloFlash.pl -lib male -CPUs 12 -read1 ../2-M4_130830_L004_R1.fastq.gz -read2 ../2-M4_130830_L004_R2.fastq.gz -dbhome /Data/filip/phyloFlash_DB/128/ -zip -readlength 100 -almosteverything -id 65 -readlimit 5000000


/opt/phyloFlash_v3.1b2/phyloFlash.pl -lib B1 -CPUs 12 -read1 ../B1_ATCACG_L002_R1_trim_001.fastq.gz -read2 ../B1_ATCACG_L002_R2_trim_001.fastq.gz -dbhome /Data/filip/phyloFlash_DB/128/ -zip -readlength 100 -almosteverything -id 65
/opt/phyloFlash_v3.1b2/phyloFlash.pl -lib B2 -CPUs 12 -read1 ../B2_CGATGT_L002_R1_trim_001.fastq.gz -read2 ../B2_CGATGT_L002_R2_trim_001.fastq.gz -dbhome /Data/filip/phyloFlash_DB/128/ -zip -readlength 100 -almosteverything -id 65
/opt/phyloFlash_v3.1b2/phyloFlash.pl -lib B4 -CPUs 12 -read1 ../B4_TTAGGC_L002_R1_trim_001.fastq.gz -read2 ../B4_TTAGGC_L002_R2_trim_001.fastq.gz -dbhome /Data/filip/phyloFlash_DB/128/ -zip -readlength 100 -almosteverything -id 65
/opt/phyloFlash_v3.1b2/phyloFlash.pl -lib B6 -CPUs 12 -read1 ../B6_TGACCA_L002_R1_trim_001.fastq.gz -read2 ../B6_TGACCA_L002_R2_trim_001.fastq.gz -dbhome /Data/filip/phyloFlash_DB/128/ -zip -readlength 100 -almosteverything -id 65
/opt/phyloFlash_v3.1b2/phyloFlash.pl -lib B7 -CPUs 12 -read1 ../B7_ACAGTG_L002_R1_trim_001.fastq.gz -read2 ../B7_ACAGTG_L002_R2_trim_001.fastq.gz -dbhome /Data/filip/phyloFlash_DB/128/ -zip -readlength 100 -almosteverything -id 65

/opt/phyloFlash_v3.1b2/phyloFlash.pl -lib G1 -CPUs 12 -read1 ../G1_CAGATC_L002_R1_trim_001.fastq.gz -read2 ../G1_CAGATC_L002_R2_trim_001.fastq.gz -dbhome /Data/filip/phyloFlash_DB/128/ -zip -readlength 100 -almosteverything -id 65
/opt/phyloFlash_v3.1b2/phyloFlash.pl -lib G2 -CPUs 12 -read1 ../G2_ACTTGA_L002_R1_trim_001.fastq.gz -read2 ../G2_ACTTGA_L002_R2_trim_001.fastq.gz -dbhome /Data/filip/phyloFlash_DB/128/ -zip -readlength 100 -almosteverything -id 65
/opt/phyloFlash_v3.1b2/phyloFlash.pl -lib G4 -CPUs 12 -read1 ../G4_GATCAG_L002_R1_trim_001.fastq.gz -read2 ../G4_GATCAG_L002_R2_trim_001.fastq.gz -dbhome /Data/filip/phyloFlash_DB/128/ -zip -readlength 100 -almosteverything -id 65
/opt/phyloFlash_v3.1b2/phyloFlash.pl -lib G6 -CPUs 12 -read1 ../G6_TAGCTT_L002_R1_trim_001.fastq.gz -read2 ../G6_TAGCTT_L002_R2_trim_001.fastq.gz -dbhome /Data/filip/phyloFlash_DB/128/ -zip -readlength 100 -almosteverything -id 65
/opt/phyloFlash_v3.1b2/phyloFlash.pl -lib G7 -CPUs 12 -read1 ../G7_GGCTAC_L002_R1_trim_001.fastq.gz -read2 ../G7_GGCTAC_L002_R2_trim_001.fastq.gz -dbhome /Data/filip/phyloFlash_DB/128/ -zip -readlength 100 -almosteverything -id 65

#DNA-Seq
/opt/phyloFlash_v3.1b2/phyloFlash.pl -lib S5_DNA -CPUs 12 -read1 ../../DNA-seq_test/s_5_1_forward_trimmed.fq.gz -read2 ../DNA-seq_test/s_5_1_reverse_trimmed.fq.gz -dbhome /Data/filip/phyloFlash_DB/128/ -zip -readlength 100 -almosteverything -id 65
/opt/phyloFlash_v3.1b2/phyloFlash.pl -lib S6_DNA -CPUs 12 -read1 ../../DNA-seq_test/s_6_1_forward_trimmed.fq.gz -read2 ../DNA-seq_test/s_6_1_reverse_trimmed.fq.gz -dbhome /Data/filip/phyloFlash_DB/128/ -zip -readlength 100 -almosteverything -id 65
```

*Transcriptome assembly in Trinity*

*Protein prediction in Transdecoder*
```
/opt/TransDecoder-TransDecoder-v5.1.0/TransDecoder.LongOrfs -t Trinity.fasta
/opt/TransDecoder-TransDecoder-v5.1.0/TransDecoder.Predict -t Trinity.fasta
```
*Species composition assesmennt in Blobtools*

*Transcript abundance estimation in RSEM*

*Differential expression analysis in EdgeR*

*Functional annotation in Trinotate*

# Arsenophonus melophagi

*Read mapping in Bowtie 2*

*Abundance estimation and transcriptome analysis in Rockhopper*

*Pseudogene annotation in Pseudo-finder*

