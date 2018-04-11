# Commands to reproduce the RNA-Seq analysis of M. ovinus published in Husnik et al.

# Melophagus ovinus
*Data quality assesment in RseQC and FastQC*

*Adapter and quality trimming in Cutadapt and Sickle*

*SSU contamination assesment in PhyloFlash*

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

