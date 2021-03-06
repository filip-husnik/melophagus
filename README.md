# Commands to reproduce RNA-Seq analyses of *Melophagus ovinus*

**Husnik, Hypsa, and Darby 2018; https://www.biorxiv.org/content/10.1101/572495v1.**

Data generated:

(1) 1 strand-specific RNA-Seq library from *M. ovinus* female

(2) 1 strand-specific RNA-Seq library from *M. ovinus* male

(3) 5 strand-specific RNA-Seq libraries for the bacteriome section of midgut (5 biological replicates)

(4) 5 strand-specific RNA-Seq libraries for the remaining sections of gut (5 biological replicates)

**Table of contents**
====================
  * [**Melophagus ovinus**](#melophagus-ovinus)
                  * [Data quality assesment in FastQC](#data-quality-assesment-in-fastqc)
                  * [Adapter and quality trimming in Cutadapt and Sickle](#adapter-and-quality-trimming-in-cutadapt-and-sickle)
                  * [SSU contamination assesment in PhyloFlash](#ssu-contamination-assesment-in-phyloflash)
                  * [Total metatranscriptome assembly in Trinity (singletons not used)](#total-metatranscriptome-assembly-in-trinity-singletons-not-used)
                  * [Completeness assesment in BUSCO (total metatranscriptome)](#completeness-assesment-in-busco-total-metatranscriptome)
                  * [Protein prediction in Transdecoder (total metatranscriptome)](#protein-prediction-in-transdecoder-total-metatranscriptome)
                  * [Species composition assesmennt for gut vs bacteriome in Blobtools](#species-composition-assesmennt-for-gut-vs-bacteriome-in-blobtools)
                  * [Transcript abundance estimation in RSEM](#transcript-abundance-estimation-in-rsem)
                  * [Relatively strict removal of contamination and lowly expressed transcripts](#relatively-strict-removal-of-contamination-and-lowly-expressed-transcripts)
                  * [Completeness assesment in BUSCO (filtered M. ovinus transcriptome)](#completeness-assesment-in-busco-filtered-m-ovinus-transcriptome)
                  * [Protein prediction in Transdecoder (filtered M. ovinus transcriptome)](#protein-prediction-in-transdecoder-filtered-m-ovinus-transcriptome)
                  * [QC samples and biological replicates](#qc-samples-and-biological-replicates)
                  * [Differential expression analysis in EdgeR](#differential-expression-analysis-in-edger)
                  * [Functional annotation in Trinotate](#functional-annotation-in-trinotate)
                  * [Create a Trinotate Web database](#create-a-trinotate-web-database)
   * [**Arsenophonus melophagi**](#arsenophonus-melophagi)
                  * [Read mapping in Bowtie 2](#read-mapping-in-bowtie-2)
                  * [Pseudogene annotation in Pseudo-finder](#pseudogene-annotation-in-pseudo-finder)
                  * [Abundance estimation and transcriptome analysis in Rockhopper](#abundance-estimation-and-transcriptome-analysis-in-rockhopper)


## Melophagus ovinus

This analysis was not run as a script, but if you decide to do so (and fix file names, locations, etc.), be cautious and use *set -e* and *set -o pipefail* in bash.

###### Data quality assesment in FastQC
```
/opt/FastQC/fastqc *.fastq.gz
```
###### Adapter and quality trimming in Cutadapt and Sickle
The data were provided preprocessed by the sequencing centre. Raw Fastq files were trimmed for the presence of Illumina adapter sequences using Cutadapt v1.1 [http://code.google.com/p/cutadapt]. Option -O 3 was used, so that the 3' end of any read which matched the adapter sequence for 3 bp or more was trimmed. The reads were further quality-trimmed by Sickle v1.200 [https://github.com/najoshi/sickle] with a minimum window quality score of 20. Reads shorter than 10 bp after trimming and singlet reads were removed.
###### SSU contamination assesment in PhyloFlash
```
#RNA-Seq
#Whole male/female libraries limited to 5 million reads
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

###### Total metatranscriptome assembly in Trinity (singletons not used)
```
/opt/trinityrnaseq-Trinity-v2.4.0/Trinity --normalize_max_read_cov 50 --full_cleanup --trimmomatic --SS_lib_type RF --CPU 6 --seqType fq --min_contig_length 300 --max_memory 150G --left B1_ATCACG_L002_R1_trim_001.fastq.gz,B2_CGATGT_L002_R1_trim_001.fastq.gz,B4_TTAGGC_L002_R1_trim_001.fastq.gz,B6_TGACCA_L002_R1_trim_001.fastq.gz,B7_ACAGTG_L002_R1_trim_001.fastq.gz,G1_CAGATC_L002_R1_trim_001.fastq.gz,G2_ACTTGA_L002_R1_trim_001.fastq.gz,G4_GATCAG_L002_R1_trim_001.fastq.gz,G6_TAGCTT_L002_R1_trim_001.fastq.gz,G7_GGCTAC_L002_R1_trim_001.fastq.gz,1-MO3_130830_L004_R1.fastq.gz,2-M4_130830_L004_R1.fastq.gz --right B1_ATCACG_L002_R2_trim_001.fastq.gz,B2_CGATGT_L002_R2_trim_001.fastq.gz,B4_TTAGGC_L002_R2_trim_001.fastq.gz,B6_TGACCA_L002_R2_trim_001.fastq.gz,B7_ACAGTG_L002_R2_trim_001.fastq.gz,G1_CAGATC_L002_R2_trim_001.fastq.gz,G2_ACTTGA_L002_R2_trim_001.fastq.gz,G4_GATCAG_L002_R2_trim_001.fastq.gz,G6_TAGCTT_L002_R2_trim_001.fastq.gz,G7_GGCTAC_L002_R2_trim_001.fastq.gz,1-MO3_130830_L004_R2.fastq.gz,2-M4_130830_L004_R2.fastq.gz
# To print basic statistics
/opt/trinityrnaseq-Trinity-v2.4.0/util/TrinityStats.pl Trinity.fasta
```
###### Completeness assesment in BUSCO (total metatranscriptome)
```
/opt/BUSCO_v3.0.2b/scripts/run_BUSCO.py -i Trinity.fasta -c 8 -o BUSCO_euk -m transcriptome -l /opt/BUSCO_v3.0.2b/databases/eukaryota_odb9/ --long
```
###### Protein prediction in Transdecoder (total metatranscriptome)
```
/opt/TransDecoder-TransDecoder-v5.1.0/TransDecoder.LongOrfs -t Trinity.fasta
/opt/TransDecoder-TransDecoder-v5.1.0/TransDecoder.Predict -t Trinity.fasta
```
###### Species composition assesmennt for gut vs bacteriome in Blobtools
```
bowtie2-build Trinity.fasta Trinity.fasta
bowtie2 -p 16 -q --mm -x Trinity.fasta -1 B1_ATCACG_L002_R1_trim_001.fastq.gz,B2_CGATGT_L002_R1_trim_001.fastq.gz,B4_TTAGGC_L002_R1_trim_001.fastq.gz,B6_TGACCA_L002_R1_trim_001.fastq.gz,B7_ACAGTG_L002_R1_trim_001.fastq.gz -2 B1_ATCACG_L002_R2_trim_001.fastq.gz,B2_CGATGT_L002_R2_trim_001.fastq.gz,B4_TTAGGC_L002_R2_trim_001.fastq.gz,B6_TGACCA_L002_R2_trim_001.fastq.gz,B7_ACAGTG_L002_R2_trim_001.fastq.gz > melophagus_bacteriome_aligned.sam

/opt/blobtools/blobtools map2cov -i Trinity.fasta -s melophagus_bacteriome_aligned.sam
rm melophagus_bacteriome_aligned.sam

bowtie2 -p 16 -q --mm -x Trinity.fasta -1 G1_CAGATC_L002_R1_trim_001.fastq.gz,G2_ACTTGA_L002_R1_trim_001.fastq.gz,G4_GATCAG_L002_R1_trim_001.fastq.gz,G6_TAGCTT_L002_R1_trim_001.fastq.gz,G7_GGCTAC_L002_R1_trim_001.fastq.gz -2 G1_CAGATC_L002_R2_trim_001.fastq.gz,G2_ACTTGA_L002_R2_trim_001.fastq.gz,G4_GATCAG_L002_R2_trim_001.fastq.gz,G6_TAGCTT_L002_R2_trim_001.fastq.gz,G7_GGCTAC_L002_R2_trim_001.fastq.gz > melophagus_gut_aligned.sam

/opt/blobtools/blobtools map2cov -i Trinity.fasta -s melophagus_gut_aligned.sam
rm melophagus_bacteriome_aligned.sam melophagus_gut_aligned.sam

bowtie2 -p 16 -q --mm -x Trinity.fasta -1 1-MO3_130830_L004_R1.fastq.gz -2 1-MO3_130830_L004_R2.fastq.gz > melophagus_female_aligned.sam
/opt/blobtools/blobtools map2cov -i Trinity.fasta -s melophagus_female_aligned.sam
rm melophagus_female_aligned.sam

bowtie2 -p 16 -q --mm -x Trinity.fasta -1 2-M4_130830_L004_R1.fastq.gz -2 2-M4_130830_L004_R2.fastq.gz > melophagus_male_aligned.sam
/opt/blobtools/blobtools map2cov -i Trinity.fasta -s melophagus_male_aligned.sam
rm melophagus_male_aligned.sam

blastn -task megablast -query Trinity.fasta -db /scratch/NCBI_NT/nt -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' -culling_limit 5 -num_threads 16 -evalue 1e-25 -max_target_seqs 5 > Trinity.fasta_assembly_vs_nt.blastn
/opt/bin/diamond blastx --query Trinity.fasta --max-target-seqs 1 --sensitive --threads 16 --db /scratch/uniprot/uniprot_ref_proteomes.diamond.dmnd --evalue 1e-25 --outfmt 6 --out Trinity.fasta.vs.uniprot_ref.mts1.1e25.out
/opt/blobtools/blobtools taxify -f Trinity.fasta.vs.uniprot_ref.mts1.1e25.out -m /scratch/uniprot/uniprot_ref_proteomes.taxids -s 0 -t 2
/opt/blobtools/blobtools create -i Trinity.fasta -t Trinity.fasta_assembly_vs_nt.blastn -t Trinity.fasta.vs.uniprot_ref.mts1.1e25.taxified.out -c melophagus_bacteriome_aligned.sam.cov -c melophagus_gut_aligned.sam.cov -c melophagus_female_aligned.sam.cov -c melophagus_male_aligned.sam.cov
/opt/blobtools/blobtools plot -i blobDB.json --rank superkingdom
/opt/blobtools/blobtools plot -i blobDB.json
/opt/blobtools/blobtools plot -i blobDB.json --format svg --noblobs
/opt/blobtools/blobtools view -i blobDB.json --rank all --hits
```

###### Transcript abundance estimation in RSEM
```
/opt/trinityrnaseq-Trinity-v2.4.0/util/align_and_estimate_abundance.pl --transcripts Trinity.fasta --seqType fq --samples_file samples_files.tsv --est_method RSEM --output_dir RSEM_abundance --aln_method bowtie2 --SS_lib_type RF --thread_count 24 --trinity_mode --prep_reference
```
###### Relatively strict removal of contamination and lowly expressed transcripts
```
/opt/trinityrnaseq-Trinity-v2.4.0/util/filter_low_expr_transcripts.pl --matrix RSEM_results/matrix.TPM.not_cross_norm --transcripts Trinity.fasta --min_expr_any 1 > Trinity_expressed_1.fasta

grep "#" -v blobDB.bestsum.table.txt | grep "Bacteria" -v | grep "Chordata" -v | grep "Kinetoplastida" -v | cut -f1 > decontaminate_transcriptome_ids.txt
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' decontaminate_transcriptome_ids.txt Trinity_expressed_1.fasta > Trinity_expressed_decontaminated.fasta
```
###### Completeness assesment in BUSCO (filtered M. ovinus transcriptome)
```
/opt/BUSCO_v3.0.2b/scripts/run_BUSCO.py -i Trinity_expressed_decontaminated.fasta -c 8 -o BUSCO_euk_M_ovinus -m transcriptome -l /opt/BUSCO_v3.0.2b/databases/eukaryota_odb9/ --long
```
###### Protein prediction in Transdecoder (filtered M. ovinus transcriptome)
```
/opt/TransDecoder-TransDecoder-v5.1.0/TransDecoder.LongOrfs -t Trinity_expressed_decontaminated.fasta
/opt/TransDecoder-TransDecoder-v5.1.0/TransDecoder.Predict -t Trinity_expressed_decontaminated.fasta
```

###### QC samples and biological replicates
```
/opt/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/PtR --matrix RSEM_results/matrix.counts.matrix --samples samples_files.tsv --CPM --log2 --min_rowSums 10 --compare_replicates
/opt/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/PtR --matrix RSEM_results/matrix.counts.matrix --min_rowSums 10 -s samples_files.tsv --log2 --CPM --sample_cor_matrix
/opt/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/PtR --matrix RSEM_results/matrix.counts.matrix -s samples_files.tsv --min_rowSums 10 --log2 --CPM --center_rows --prin_comp 3
```

###### Differential expression analysis in EdgeR
```
/opt/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix RSEM_results/matrix.counts.matrix --method edgeR --samples_file samples_files.tsv
```

###### Functional annotation in Trinotate
```
blastx -query Trinity.fasta -db /opt/Trinotate-Trinotate-v3.1.1/uniprot_sprot.pep -num_threads 16 -max_target_seqs 1 -outfmt 6 > blastx.outfmt6
blastp -query Trinity.fasta.transdecoder.pep -db /opt/Trinotate-Trinotate-v3.1.1/uniprot_sprot.pep -num_threads 16 -max_target_seqs 1 -outfmt 6 > blastp.outfmt6
hmmscan --cpu 16 --domtblout TrinotatePFAM.out /opt/Trinotate-Trinotate-v3.1.1/Pfam-A.hmm Trinity.fasta.transdecoder.pep > pfam.log
/opt/signalp-4.1/signalp_4.1 -f short -n signalp.out Trinity.fasta.transdecoder.pep
/opt/tmhmm-2.0c/bin/tmhmm --short < Trinity.fasta.transdecoder.pep > tmhmm.out

# Only eukaryotic rRNAs
#/opt/Trinotate-Trinotate-v3.1.1/util/rnammer_support/RnammerTranscriptome.pl --transcriptome Trinity.fasta --path_to_rnammer /opt/rnammer-1.2/rnammer
mv Trinity.fasta.rnammer.gff Trinity.euk.fasta.rnammer.gff

# To also get bacterial rRNA coordinates
#perl /opt/rnammer-1.2/rnammer -S bac -m tsu,lsu,ssu -gff bac.tmp.superscaff.rnammer.gff < transcriptSuperScaffold.fasta
#/opt/Trinotate-Trinotate-v3.1.1/util/rnammer_support/util/rnammer_supperscaffold_gff_to_indiv_transcripts.pl -R bac.tmp.superscaff.rnammer.gff -T transcriptSuperScaffold.bed > Trinity.bac.fasta.rnammer.gff

#/opt/trinityrnaseq-Trinity-v2.4.0/util/support_scripts/get_Trinity_gene_to_trans_map.pl Trinity.fasta >  Trinity.fasta.gene_trans_map

/opt/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite init --gene_trans_map Trinity.fasta.gene_trans_map --transcript_fasta Trinity.fasta --transdecoder_pep Trinity.fasta.transdecoder.pep
/opt/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite LOAD_swissprot_blastp blastp.outfmt6
/opt/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite LOAD_swissprot_blastx blastx.outfmt6
/opt/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite LOAD_pfam TrinotatePFAM.out
/opt/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite LOAD_tmhmm tmhmm.out
/opt/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite LOAD_signalp signalp.out
/opt/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite LOAD_rnammer Trinity.euk.fasta.rnammer.gff
/opt/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite LOAD_rnammer Trinity.bac.fasta.rnammer.gff
/opt/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite report > trinotate_annotation_report.xls
```

###### Create a Trinotate Web database

```
# First cp the Trinotate.sqlite database into a new folder and name it TrinotateWeb.sqlite (it's better to keep it separated)
cp Trinotate.sqlite TrinotateWeb/TrinotateWeb.sqlite
cd TrinotateWeb/

# Import the fpkm and DE analysis stuff (matrix.counts.matrix.cond_A_vs_cond_B.edgeR.DE_results)
/opt/Trinotate-Trinotate-v3.1.1/util/transcript_expression/import_expression_and_DE_results.pl \
          --sqlite TrinotateWeb.sqlite \
          --samples_file samples_files.tsv \
          --count_matrix matrix.counts.matrix \
          --fpkm_matrix matrix.TMM.EXPR.matrix \
          --DE_dir edgeR_DE_results/ \
          --transcript_mode
          
# Import text annotations
/opt/Trinotate-Trinotate-v3.1.1/util/annotation_importer/import_transcript_names.pl TrinotateWeb.sqlite trinotate_annotation_report.xls

# Alternatively, import species assignments from the Blobtools taxonomy table
#gene_id (tab) transcript_id (tab) annotation text
#create a taxonomy file first
#grep -v "#" blobDB.bestsum.table.txt | cut -f 8,12,16,20,24,28 | sed s"/\t/:/"g > species_annotations.txt
#grep -v "#" blobDB.bestsum.table.txt | cut -f 1 > transcripts.txt
#cut -f 1,2,3,4 transcripts.txt -d '_' > genes.txt
#paste genes.txt transcripts.txt species_annotations.txt > taxonomy_annotation.tsv
#rm species_annotations.txt transcripts.txt genes.txt

#/opt/Trinotate-Trinotate-v3.1.1/util/annotation_importer/import_transcript_names.pl TrinotateWeb.sqlite taxonomy_annotation.tsv

#Run the webserver on your own computer
#/home/filip/programs/Trinotate-Trinotate-v3.1.1/run_TrinotateWebserver.pl 8080
#http://localhost:8080/cgi-bin/index.cgi
```

## Arsenophonus melophagi

###### Read mapping in Bowtie 2
```
bowtie2-build ARM.fasta ARM.fasta.
bowtie2 -a --no-unal --threads 8 -x ARM.fasta -1 B1_ATCACG_L002_R1_trim_001.fastq,B2_CGATGT_L002_R1_trim_0 01.fastq,B4_TTAGGC_L002_R1_trim_001.fastq,B6_TGACCA_L002_R1_trim_001.fastq,B7_ACAGTG_L002_R1_trim_001.fastq -2 B1_ATCACG_L002_R2_trim_001.fastq,B2_CGATGT_L002_R2_trim_001.fastq,B4_TTAGGC_L002_R2_trim_001.fastq,B6_TGACCA_L002_R2_trim_001.fastq,B7_ACAGTG_L002_R2_trim_001.fastq  | samtools view -bS - > ARM_PE.bam
samtools sort ARM_PE.bam ARM_PE_sorted
samtools index ARM_PE_sorted.bam ARM_PE_sorted.bam.bai
```

###### Pseudogene annotation in Pseudo-finder

Initial annotation in PROKKA
```
prokka --compliant --rfam --gram neg AME.fasta
```

Pseudo-finder [https://github.com/filip-husnik/pseudo-finder] pseudogene annotation/visualization and manual curation in Artemis/Bamview according to expression data.

```
python3 pseudo_finder.py annotate --genome AME.gbf --outprefix AMEps --database nr --outformat 124 --threads 18
python3 pseudo_finder.py visualize --genome AME.gbf --outprefix AMEps_plots --blastp AMEps_AME.gbf_proteome.faa.blastP_output.tsv --blastx AMEps_AME.gbf_intergenic.fasta.blastX_output.tsv --hitcap 15
python3 pseudo_finder.py map --genome AME.gbf --gff AMEps_AME.gbf_pseudos.gff --outprefix AMEps_map
```

###### Abundance estimation and transcriptome analysis in Rockhopper
Rockhopper has a graphical user interface and its default analysis with five biological replicates (bacteriome libraries) was used for *Arsenophonus melophagi* [https://cs.wellesley.edu/~btjaden/Rockhopper/].

