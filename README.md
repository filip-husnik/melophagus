# Commands to reproduce the RNA-Seq analysis of M. ovinus published in Husnik et al.

# Melophagus ovinus

```
#!/usr/bin/bash
set -e
set -o pipefail
```

*Data quality assesment in FastQC*
```
/opt/FastQC/fastqc *.fastq.gz
```
*Adapter and quality trimming in Cutadapt and Sickle*
```
# Carried out by the sequencing centre.
```
*SSU contamination assesment in PhyloFlash*
```
#RNA-Seq
/opt/phyloFlash_v3.1b2/phyloFlash.pl -lib female -CPUs 12 -read1 ../1-MO3_130830_L004_R1.fastq.gz -read2 ../1-MO3_130830_L004_R2.fastq.gz -dbhome /Data/filip/phyloFlash_DB/128/ -zip -readlength 100 -almosteverything -id 65 -readlimit 5000000
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

*Total metatranscriptome assembly in Trinity (singletons not used)*
```
/opt/trinityrnaseq-Trinity-v2.4.0/Trinity --normalize_max_read_cov 50 --full_cleanup --trimmomatic --SS_lib_type RF --CPU 6 --seqType fq --min_contig_length 300 --max_memory 150G --left B1_ATCACG_L002_R1_trim_001.fastq.gz,B2_CGATGT_L002_R1_trim_001.fastq.gz,B4_TTAGGC_L002_R1_trim_001.fastq.gz,B6_TGACCA_L002_R1_trim_001.fastq.gz,B7_ACAGTG_L002_R1_trim_001.fastq.gz,G1_CAGATC_L002_R1_trim_001.fastq.gz,G2_ACTTGA_L002_R1_trim_001.fastq.gz,G4_GATCAG_L002_R1_trim_001.fastq.gz,G6_TAGCTT_L002_R1_trim_001.fastq.gz,G7_GGCTAC_L002_R1_trim_001.fastq.gz,1-MO3_130830_L004_R1.fastq.gz,2-M4_130830_L004_R1.fastq.gz --right B1_ATCACG_L002_R2_trim_001.fastq.gz,B2_CGATGT_L002_R2_trim_001.fastq.gz,B4_TTAGGC_L002_R2_trim_001.fastq.gz,B6_TGACCA_L002_R2_trim_001.fastq.gz,B7_ACAGTG_L002_R2_trim_001.fastq.gz,G1_CAGATC_L002_R2_trim_001.fastq.gz,G2_ACTTGA_L002_R2_trim_001.fastq.gz,G4_GATCAG_L002_R2_trim_001.fastq.gz,G6_TAGCTT_L002_R2_trim_001.fastq.gz,G7_GGCTAC_L002_R2_trim_001.fastq.gz,1-MO3_130830_L004_R2.fastq.gz,2-M4_130830_L004_R2.fastq.gz
# To print basic statistics
/opt/trinityrnaseq-Trinity-v2.4.0/util/TrinityStats.pl Trinity.fasta
```
*Completeness assesment in BUSCO (total metatranscriptome)*
```
/opt/BUSCO_v3.0.2b/scripts/run_BUSCO.py -i Trinity.fasta -c 8 -o BUSCO_euk -m transcriptome -l /opt/BUSCO_v3.0.2b/databases/eukaryota_odb9/ --long
```
*Protein prediction in Transdecoder (total metatranscriptome)*
```
/opt/TransDecoder-TransDecoder-v5.1.0/TransDecoder.LongOrfs -t Trinity.fasta
/opt/TransDecoder-TransDecoder-v5.1.0/TransDecoder.Predict -t Trinity.fasta
```
*Species composition assesmennt for gut vs bacteriome in Blobtools*
```
bowtie2-build Trinity.fasta Trinity.fasta
bowtie2 -p 16 -q --mm -x Trinity.fasta -1 B1_ATCACG_L002_R1_trim_001.fastq.gz,B2_CGATGT_L002_R1_trim_001.fastq.gz,B4_TTAGGC_L002_R1_trim_001.fastq.gz,B6_TGACCA_L002_R1_trim_001.fastq.gz,B7_ACAGTG_L002_R1_trim_001.fastq.gz -2 B1_ATCACG_L002_R2_trim_001.fastq.gz,B2_CGATGT_L002_R2_trim_001.fastq.gz,B4_TTAGGC_L002_R2_trim_001.fastq.gz,B6_TGACCA_L002_R2_trim_001.fastq.gz,B7_ACAGTG_L002_R2_trim_001.fastq.gz > melophagus_bacteriome_aligned.sam

/opt/blobtools/blobtools map2cov -i Trinity.fasta -s melophagus_bacteriome_aligned.sam
rm melophagus_bacteriome_aligned.sam

bowtie2 -p 16 -q --mm -x Trinity.fasta -1 G1_CAGATC_L002_R1_trim_001.fastq.gz,G2_ACTTGA_L002_R1_trim_001.fastq.gz,G4_GATCAG_L002_R1_trim_001.fastq.gz,G6_TAGCTT_L002_R1_trim_001.fastq.gz,G7_GGCTAC_L002_R1_trim_001.fastq.gz -2 G1_CAGATC_L002_R2_trim_001.fastq.gz,G2_ACTTGA_L002_R2_trim_001.fastq.gz,G4_GATCAG_L002_R2_trim_001.fastq.gz,G6_TAGCTT_L002_R2_trim_001.fastq.gz,G7_GGCTAC_L002_R2_trim_001.fastq.gz > melophagus_gut_aligned.sam

/opt/blobtools/blobtools map2cov -i Trinity.fasta -s melophagus_gut_aligned.sam
rm melophagus_bacteriome_aligned.sam melophagus_gut_aligned.sam

blastn -task megablast -query Trinity.fasta -db /scratch/NCBI_NT/nt -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' -culling_limit 5 -num_threads 16 -evalue 1e-25 -max_target_seqs 5 > Trinity.fasta_assembly_vs_nt.blastn
/opt/bin/diamond blastx --query Trinity.fasta --max-target-seqs 1 --sensitive --threads 16 --db /scratch/uniprot/uniprot_ref_proteomes.diamond.dmnd --evalue 1e-25 --outfmt 6 --out Trinity.fasta.vs.uniprot_ref.mts1.1e25.out
/opt/blobtools/blobtools taxify -f Trinity.fasta.vs.uniprot_ref.mts1.1e25.out -m /scratch/uniprot/uniprot_ref_proteomes.taxids -s 0 -t 2
/opt/blobtools/blobtools create -i Trinity.fasta -t Trinity.fasta_assembly_vs_nt.blastn -t Trinity.fasta.vs.uniprot_ref.mts1.1e25.taxified.out -c melophagus_bacteriome_aligned.sam.cov -c melophagus_gut_aligned.sam.cov
/opt/blobtools/blobtools plot -i blobDB.json --rank superkingdom
/opt/blobtools/blobtools plot -i blobDB.json
/opt/blobtools/blobtools view -i blobDB.json --rank all --hits
```

*Transcript abundance estimation in RSEM*
```
/opt/trinityrnaseq-Trinity-v2.4.0/util/align_and_estimate_abundance.pl --transcripts Trinity.fasta --seqType fq --samples_file samples_files.tsv --est_method RSEM --output_dir RSEM_abundance --aln_method bowtie2 --SS_lib_type RF --thread_count 24 --trinity_mode --prep_reference
```
*Relatively strict removal of contamination and lowly expressed transcripts*
```
/opt/trinityrnaseq-Trinity-v2.4.0/util/filter_low_expr_transcripts.pl --matrix RSEM_results/matrix.TPM.not_cross_norm --transcripts Trinity.fasta --min_expr_any 1 > Trinity_expressed_1.fasta

grep "#" -v blobDB.bestsum.table.txt | grep "Bacteria" -v | grep "Chordata" -v | grep "Kinetoplastida" -v | cut -f1 > decontaminate_transcriptome_ids.txt
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' decontaminate_transcriptome_ids.txt Trinity_expressed_1.fasta > Trinity_expressed_decontaminated.fasta
```
*Completeness assesment in BUSCO (filtered M. ovinus transcriptome)*
```
/opt/BUSCO_v3.0.2b/scripts/run_BUSCO.py -i Trinity_expressed_decontaminated.fasta -c 8 -o BUSCO_euk_M_ovinus -m transcriptome -l /opt/BUSCO_v3.0.2b/databases/eukaryota_odb9/ --long
```
*Protein prediction in Transdecoder (filtered M. ovinus transcriptome)*
```
/opt/TransDecoder-TransDecoder-v5.1.0/TransDecoder.LongOrfs -t Trinity_expressed_decontaminated.fasta
/opt/TransDecoder-TransDecoder-v5.1.0/TransDecoder.Predict -t Trinity_expressed_decontaminated.fasta
```


*QC samples and biological replicates*
```
/opt/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/PtR --matrix RSEM_results/matrix.counts.matrix --samples samples_files.tsv --CPM --log2 --min_rowSums 10 --compare_replicates

/opt/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/PtR --matrix RSEM_results/matrix.counts.matrix --min_rowSums 10 -s samples_files.tsv --log2 --CPM --sample_cor_matrix

/opt/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/PtR --matrix RSEM_results/matrix.counts.matrix -s samples_files.tsv --min_rowSums 10 --log2 --CPM --center_rows --prin_comp 3
```

*Differential expression analysis in EdgeR*
```
/opt/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix RSEM_results/matrix.counts.matrix --method edgeR --samples_file samples_files.tsv
```

*Functional annotation in Trinotate*
```
blastx -query Trinity.fasta -db /opt/Trinotate-Trinotate-v3.1.1/uniprot_sprot.pep -num_threads 16 -max_target_seqs 1 -outfmt 6 > blastx.outfmt6
blastp -query Trinity.fasta.transdecoder.pep -db /opt/Trinotate-Trinotate-v3.1.1/uniprot_sprot.pep -num_threads 16 -max_target_seqs 1 -outfmt 6 > blastp.outfmt6
hmmscan --cpu 16 --domtblout TrinotatePFAM.out /opt/Trinotate-Trinotate-v3.1.1/Pfam-A.hmm Trinity.fasta.transdecoder.pep > pfam.log
/opt/signalp-4.1/signalp_4.1 -f short -n signalp.out Trinity.fasta.transdecoder.pep
/opt/tmhmm-2.0c/bin/tmhmm --short < Trinity.fasta.transdecoder.pep > tmhmm.out

# Only eukaryotic rRNAs (on Rosetta)
#/opt/Trinotate-Trinotate-v3.1.1/util/rnammer_support/RnammerTranscriptome.pl --transcriptome Trinity.fasta --path_to_rnammer /opt/rnammer-1.2/rnammer 

# To also get bacterial rRNA coordinates (on Rosetta)
#perl /opt/rnammer-1.2/rnammer -S bac -m tsu,lsu,ssu -gff bac.tmp.superscaff.rnammer.gff < transcriptSuperScaffold.fasta
#/opt/Trinotate-Trinotate-v3.1.1/util/rnammer_support/util/rnammer_supperscaffold_gff_to_indiv_transcripts.pl -R bac.tmp.superscaff.rnammer.gff -T transcriptSuperScaffold.bed > Trinity.bac.fasta.rnammer.gff

#/opt/trinityrnaseq-Trinity-v2.4.0/util/support_scripts/get_Trinity_gene_to_trans_map.pl Trinity.fasta >  Trinity.fasta.gene_trans_map

/opt/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite init --gene_trans_map Trinity.fasta.gene_trans_map --transcript_fasta Trinity.fasta --transdecoder_pep Trinity.fasta.transdecoder.pep
/opt/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite LOAD_swissprot_blastp blastp.outfmt6
/opt/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite LOAD_swissprot_blastx blastx.outfmt6
/opt/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite LOAD_pfam TrinotatePFAM.out
/opt/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite LOAD_tmhmm tmhmm.out
/opt/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite LOAD_signalp signalp.out
/opt/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite LOAD_rnammer Trinity.fasta.rnammer.gff
/opt/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite report > trinotate_annotation_report.xls
```

*Create a Trinotate Web database (not properly tested)*

```
# First create a boiler plate database TrinotateWeb.sqlite

#MISSING!

# Import the fpkm and DE analysis stuff
/opt/Trinotate-Trinotate-v3.1.1/util/transcript_expression/import_expression_and_DE_results.pl \
          --sqlite TrinotateWeb.sqlite \
          --samples_file samples_n_reads_described.txt \
          --count_matrix Trinity_trans.counts.matrix \
          --fpkm_matrix Trinity_trans.counts.matrix.TMM_normalized.FPKM \
          --DE_dir edgeR_trans/ \
          --transcript_mode
          
# Import text annotations
/opt/Trinotate-Trinotate-v3.1.1/util/annotation_importer/import_transcript_names.pl TrinotateWeb.sqlite trinotate_annotation_report.xls

# Import species assignment from the Blobtools taxonomy table
#gene_id (tab) transcript_id (tab) annotation text

#MISSING!

#Run the webserver
/opt/Trinotate-Trinotate-v3.1.1/run_TrinotateWebserver.pl 8080
#http://localhost:8080/cgi-bin/index.cgi
```

# Arsenophonus melophagi

*Read mapping in Bowtie 2*

*Abundance estimation and transcriptome analysis in Rockhopper*

*Pseudogene annotation in Pseudo-finder*

