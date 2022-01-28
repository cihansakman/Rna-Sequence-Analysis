# Rna-Sequence-Analysis


## Loading and Overview of Data
```
prefetch SRR
fastq-dump SRR.sra --split-files
fastqc *.fastq
multiqc .
```

## Adapter Trim
```
$RNA_HOME/student_tools/flexbar-3.4.0-linux/flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_REFS_DIR/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_MY_DATA/SRR3474918/SRR3474918_1.fastq --reads2 $RNA_MY_DATA/SRR3474918/SRR3474918_2.fastq --target $RNA_MY_DATA_TRIM/SRR3474918	
$RNA_HOME/student_tools/flexbar-3.4.0-linux/flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_REFS_DIR/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_MY_DATA/SRR3475345/SRR3475345_1.fastq --reads2 $RNA_MY_DATA/SRR3475345/SRR3475345_2.fastq --target $RNA_MY_DATA_TRIM/SRR3475345	
$RNA_HOME/student_tools/flexbar-3.4.0-linux/flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_REFS_DIR/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_MY_DATA/SRR3474779/SRR3474779_1.fastq --reads2 $RNA_MY_DATA/SRR3474779/SRR3474779_2.fastq --target $RNA_MY_DATA_TRIM/SRR3474779	
$RNA_HOME/student_tools/flexbar-3.4.0-linux/flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_REFS_DIR/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_MY_DATA/SRR3475326/SRR3475326_1.fastq --reads2 $RNA_MY_DATA/SRR3475326/SRR3475326_2.fastq --target $RNA_MY_DATA_TRIM/SRR3475326	
$RNA_HOME/student_tools/flexbar-3.4.0-linux/flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_REFS_DIR/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_MY_DATA/SRR3474795/SRR3474795_1.fastq --reads2 $RNA_MY_DATA/SRR3474795/SRR3474795_2.fastq --target $RNA_MY_DATA_TRIM/SRR3474795	
$RNA_HOME/student_tools/flexbar-3.4.0-linux/flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_REFS_DIR/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_MY_DATA/SRR3475332/SRR3475332_1.fastq --reads2 $RNA_MY_DATA/SRR3475332/SRR3475332_2.fastq --target $RNA_MY_DATA_TRIM/SRR3475332	
$RNA_HOME/student_tools/flexbar-3.4.0-linux/flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_REFS_DIR/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_MY_DATA/SRR3474972/SRR3474972_1.fastq --reads2 $RNA_MY_DATA/SRR3474972/SRR3474972_2.fastq --target $RNA_MY_DATA_TRIM/SRR3474972	
$RNA_HOME/student_tools/flexbar-3.4.0-linux/flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_REFS_DIR/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_MY_DATA/SRR3475362/SRR3475362_1.fastq --reads2 $RNA_MY_DATA/SRR3475362/SRR3475362_2.fastq --target $RNA_MY_DATA_TRIM/SRR3475362	
$RNA_HOME/student_tools/flexbar-3.4.0-linux/flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_REFS_DIR/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_MY_DATA/SRR3475081/SRR3475081_1.fastq --reads2 $RNA_MY_DATA/SRR3475081/SRR3475081_2.fastq --target $RNA_MY_DATA_TRIM/SRR3475081	
$RNA_HOME/student_tools/flexbar-3.4.0-linux/flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_REFS_DIR/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_MY_DATA/SRR3475375/SRR3475375_1.fastq --reads2 $RNA_MY_DATA/SRR3475375/SRR3475375_2.fastq --target $RNA_MY_DATA_TRIM/SRR3475375	
$RNA_HOME/student_tools/flexbar-3.4.0-linux/flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_REFS_DIR/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_MY_DATA/SRR3475089/SRR3475089_1.fastq --reads2 $RNA_MY_DATA/SRR3475089/SRR3475089_2.fastq --target $RNA_MY_DATA_TRIM/SRR3475089	
$RNA_HOME/student_tools/flexbar-3.4.0-linux/flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_REFS_DIR/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_MY_DATA/SRR3475377/SRR3475377_1.fastq --reads2 $RNA_MY_DATA/SRR3475377/SRR3475377_2.fastq --target $RNA_MY_DATA_TRIM/SRR3475377
```

## Alignment
```
#LUSC ALIGNMENTS
hisat2 -p 16 --rg-id=SRR3474779 --rg SM:LUSC --rg LB:SRR3474779-LUSC --summary-file SRR3474779.txt --rg PL:ILLUMINA -x $RNA_MY_REF_INDEX --dta --rna-strandness RF -1 $RNA_MY_DATA_TRIM/SRR3474779_1.fastq.gz -2 $RNA_MY_DATA_TRIM/SRR3474779_2.fastq.gz -S ./SRR3474779.sam
hisat2 -p 16 --rg-id=SRR3474795 --rg SM:LUSC --rg LB:SRR3474795-LUSC --summary-file SRR3474795.txt --rg PL:ILLUMINA -x $RNA_MY_REF_INDEX --dta --rna-strandness RF -1 $RNA_MY_DATA_TRIM/SRR3474795_1.fastq.gz -2 $RNA_MY_DATA_TRIM/SRR3474795_2.fastq.gz -S ./SRR3474795.sam
hisat2 -p 16 --rg-id=SRR3474918 --rg SM:LUSC --rg LB:SRR3474918-LUSC --summary-file SRR3474918.txt --rg PL:ILLUMINA -x $RNA_MY_REF_INDEX --dta --rna-strandness RF -1 $RNA_MY_DATA_TRIM/SRR3474918_1.fastq.gz -2 $RNA_MY_DATA_TRIM/SRR3474918_2.fastq.gz -S ./SRR3474918.sam
hisat2 -p 16 --rg-id=SRR3474972 --rg SM:LUSC --rg LB:SRR3474972-LUSC --summary-file SRR3474972.txt --rg PL:ILLUMINA -x $RNA_MY_REF_INDEX --dta --rna-strandness RF -1 $RNA_MY_DATA_TRIM/SRR3474972_1.fastq.gz -2 $RNA_MY_DATA_TRIM/SRR3474972_2.fastq.gz -S ./SRR3474972.sam
hisat2 -p 16 --rg-id=SRR3475081 --rg SM:LUSC --rg LB:SRR3475081-LUSC --summary-file SRR3475081.txt --rg PL:ILLUMINA -x $RNA_MY_REF_INDEX --dta --rna-strandness RF -1 $RNA_MY_DATA_TRIM/SRR3475081_1.fastq.gz -2 $RNA_MY_DATA_TRIM/SRR3475081_2.fastq.gz -S ./SRR3475081.sam
hisat2 -p 16 --rg-id=SRR3475089 --rg SM:LUSC --rg LB:SRR3475089-LUSC --summary-file SRR3475089.txt --rg PL:ILLUMINA -x $RNA_MY_REF_INDEX --dta --rna-strandness RF -1 $RNA_MY_DATA_TRIM/SRR3475089_1.fastq.gz -2 $RNA_MY_DATA_TRIM/SRR3475089_2.fastq.gz -S ./SRR3475089.sam
#CONTROL ALIGNMENTS
hisat2 -p 16 --rg-id=SRR3475326 --rg SM:NORMAL --rg LB:SRR3475326-NORMAL --summary-file SRR3475326.txt --rg PL:ILLUMINA -x $RNA_MY_REF_INDEX --dta --rna-strandness RF -1 $RNA_MY_DATA_TRIM/SRR3475326_1.fastq.gz -2 $RNA_MY_DATA_TRIM/SRR3475326_2.fastq.gz -S ./SRR3475326.sam
hisat2 -p 16 --rg-id=SRR3475332 --rg SM:NORMAL --rg LB:SRR3475332-NORMAL --summary-file SRR3475332.txt --rg PL:ILLUMINA -x $RNA_MY_REF_INDEX --dta --rna-strandness RF -1 $RNA_MY_DATA_TRIM/SRR3475332_1.fastq.gz -2 $RNA_MY_DATA_TRIM/SRR3475332_2.fastq.gz -S ./SRR3475332.sam
hisat2 -p 16 --rg-id=SRR3475345 --rg SM:NORMAL --rg LB:SRR3475345-NORMAL --summary-file SRR3475345.txt --rg PL:ILLUMINA -x $RNA_MY_REF_INDEX --dta --rna-strandness RF -1 $RNA_MY_DATA_TRIM/SRR3475345_1.fastq.gz -2 $RNA_MY_DATA_TRIM/SRR3475345_2.fastq.gz -S ./SRR3475345.sam
hisat2 -p 16 --rg-id=SRR3475362 --rg SM:NORMAL --rg LB:SRR3475362-NORMAL --summary-file SRR3475362.txt --rg PL:ILLUMINA -x $RNA_MY_REF_INDEX --dta --rna-strandness RF -1 $RNA_MY_DATA_TRIM/SRR3475362_1.fastq.gz -2 $RNA_MY_DATA_TRIM/SRR3475362_2.fastq.gz -S ./SRR3475362.sam
hisat2 -p 16 --rg-id=SRR3475375 --rg SM:NORMAL --rg LB:SRR3475375-NORMAL --summary-file SRR3475375.txt --rg PL:ILLUMINA -x $RNA_MY_REF_INDEX --dta --rna-strandness RF -1 $RNA_MY_DATA_TRIM/SRR3475375_1.fastq.gz -2 $RNA_MY_DATA_TRIM/SRR3475375_2.fastq.gz -S ./SRR3475375.sam
hisat2 -p 16 --rg-id=SRR3475377 --rg SM:NORMAL --rg LB:SRR3475377-NORMAL --summary-file SRR3475377.txt --rg PL:ILLUMINA -x $RNA_MY_REF_INDEX --dta --rna-strandness RF -1 $RNA_MY_DATA_TRIM/SRR3475377_1.fastq.gz -2 $RNA_MY_DATA_TRIM/SRR3475377_2.fastq.gz -S ./SRR3475377.sam
```
### SAM to BAM Conversion
```
samtools sort -@ 16 -o SRR3474779.bam SRR3474779.sam
samtools sort -@ 16 -o SRR3474795.bam SRR3474795.sam
samtools sort -@ 16 -o SRR3474918.bam SRR3474918.sam
samtools sort -@ 16 -o SRR3474972.bam SRR3474972.sam
samtools sort -@ 16 -o SRR3475081.bam SRR3475081.sam
samtools sort -@ 16 -o SRR3475089.bam SRR3475089.sam
samtools sort -@ 16 -o SRR3475326.bam SRR3475326.sam
samtools sort -@ 16 -o SRR3475332.bam SRR3475332.sam
samtools sort -@ 16 -o SRR3475345.bam SRR3475345.sam
samtools sort -@ 16 -o SRR3475362.bam SRR3475362.sam
samtools sort -@ 16 -o SRR3475375.bam SRR3475375.sam
samtools sort -@ 16 -o SRR3475377.bam SRR3475377.sam
```

### Merge BAM Files
```
#################### Merge HISAT2 BAM files ####################
cd $RNA_HOME/my_alignments/hisat2
java -Xmx2g -jar $RNA_HOME/student_tools/picard.jar MergeSamFiles OUTPUT=LUSC.bam INPUT=SRR3474779.bam INPUT=SRR3474795.bam INPUT=SRR3474918.bam INPUT=SRR3474972.bam INPUT=SRR3475081.bam INPUT=SRR3475089.bam
java -Xmx2g -jar $RNA_HOME/student_tools/picard.jar MergeSamFiles OUTPUT=NORMAL.bam INPUT=SRR3475326.bam INPUT=SRR3475332.bam INPUT=SRR3475345.bam INPUT=SRR3475362.bam INPUT=SRR3475375.bam INPUT=SRR3475377.bam
```
## Expression
```
#LUSC SAMPLES
stringtie -p 16 -G $RNA_MY_REF_GTF -e -B -o SRR3474779/transcripts.gtf -A SRR3474779/gene_abundances.tsv $RNA_MY_ALIGN_DIR/SRR3474779.bam
stringtie -p 16 -G $RNA_MY_REF_GTF -e -B -o SRR3474795/transcripts.gtf -A SRR3474795/gene_abundances.tsv $RNA_MY_ALIGN_DIR/SRR3474795.bam
stringtie -p 16 -G $RNA_MY_REF_GTF -e -B -o SRR3474918/transcripts.gtf -A SRR3474918/gene_abundances.tsv $RNA_MY_ALIGN_DIR/SRR3474918.bam
stringtie -p 16 -G $RNA_MY_REF_GTF -e -B -o SRR3474972/transcripts.gtf -A SRR3474972/gene_abundances.tsv $RNA_MY_ALIGN_DIR/SRR3474972.bam
stringtie -p 16 -G $RNA_MY_REF_GTF -e -B -o SRR3475081/transcripts.gtf -A SRR3475081/gene_abundances.tsv $RNA_MY_ALIGN_DIR/SRR3475081.bam
stringtie -p 16 -G $RNA_MY_REF_GTF -e -B -o SRR3475089/transcripts.gtf -A SRR3475089/gene_abundances.tsv $RNA_MY_ALIGN_DIR/SRR3475089.bam
#NORMAL SAMPLES
stringtie -p 16 -G $RNA_MY_REF_GTF -e -B -o SRR3475326/transcripts.gtf -A SRR3475326/gene_abundances.tsv $RNA_MY_ALIGN_DIR/SRR3475326.bam
stringtie -p 16 -G $RNA_MY_REF_GTF -e -B -o SRR3475332/transcripts.gtf -A SRR3475332/gene_abundances.tsv $RNA_MY_ALIGN_DIR/SRR3475332.bam
stringtie -p 16 -G $RNA_MY_REF_GTF -e -B -o SRR3475345/transcripts.gtf -A SRR3475345/gene_abundances.tsv $RNA_MY_ALIGN_DIR/SRR3475345.bam
stringtie -p 16 -G $RNA_MY_REF_GTF -e -B -o SRR3475362/transcripts.gtf -A SRR3475362/gene_abundances.tsv $RNA_MY_ALIGN_DIR/SRR3475362.bam
stringtie -p 16 -G $RNA_MY_REF_GTF -e -B -o SRR3475375/transcripts.gtf -A SRR3475375/gene_abundances.tsv $RNA_MY_ALIGN_DIR/SRR3475375.bam
stringtie -p 16 -G $RNA_MY_REF_GTF -e -B -o SRR3475377/transcripts.gtf -A SRR3475377/gene_abundances.tsv $RNA_MY_ALIGN_DIR/SRR3475377.bam
```
In the final, we can merge all the samples in a file and observe the FPKM and TPM values with following codes.
```
wget
https://raw.githubusercontent.com/griffithlab/rnaseq_tutorial/master/scripts/stringtie_expression_matrix.pl
chmod +x stringtie_expression_matrix.pl
./stringtie_expression_matrix.pl --expression_metric=TPM --result_dirs='SRR3474779, SRR3474795, SRR3474918, SRR3474972, SRR3475081, SRR3475089, SRR3475326, SRR3475332, SRR3475345, SRR3475362, SRR3475375, SRR3475377' --transcript_matrix_file=transcript_tpm_all_samples.tsv --gene_matrix_file=gene_tpm_all_samples.tsv
./stringtie_expression_matrix.pl --expression_metric=FPKM --result_dirs='SRR3474779, SRR3474795, SRR3474918, SRR3474972, SRR3475081, SRR3475089, SRR3475326, SRR3475332, SRR3475345, SRR3475362, SRR3475375, SRR3475377' --transcript_matrix_file=transcript_fpkm_all_samples.tsv --gene_matrix_file=gene_fpkm_all_samples.tsv
./stringtie_expression_matrix.pl --expression_metric=Coverage --result_dirs='SRR3474779, SRR3474795, SRR3474918, SRR3474972, SRR3475081, SRR3475089, SRR3475326, SRR3475332, SRR3475345, SRR3475362, SRR3475375, SRR3475377' --transcript_matrix_file=transcript_coverage_all_samples.tsv --gene_matrix_file=gene_coverage_all_samples.tsv
```
## HTSEQ-Count
```
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $RNA_MY_ALIGN_DIR/SRR3474779.bam $RNA_MY_REF_GTF > SRR3474779.tsv
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $RNA_MY_ALIGN_DIR/SRR3474795.bam $RNA_MY_REF_GTF > SRR3474795.tsv
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $RNA_MY_ALIGN_DIR/SRR3474918.bam $RNA_MY_REF_GTF > SRR3474918.tsv
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $RNA_MY_ALIGN_DIR/SRR3474972.bam $RNA_MY_REF_GTF > SRR3474972.tsv
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $RNA_MY_ALIGN_DIR/SRR3475081.bam $RNA_MY_REF_GTF > SRR3475081.tsv
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $RNA_MY_ALIGN_DIR/SRR3475089.bam $RNA_MY_REF_GTF > SRR3475089.tsv
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $RNA_MY_ALIGN_DIR/SRR3475326.bam $RNA_MY_REF_GTF > SRR3475326.tsv
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $RNA_MY_ALIGN_DIR/SRR3475332.bam $RNA_MY_REF_GTF > SRR3475332.tsv
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $RNA_MY_ALIGN_DIR/SRR3475345.bam $RNA_MY_REF_GTF > SRR3475345.tsv
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $RNA_MY_ALIGN_DIR/SRR3475362.bam $RNA_MY_REF_GTF > SRR3475362.tsv
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $RNA_MY_ALIGN_DIR/SRR3475375.bam $RNA_MY_REF_GTF > SRR3475375.tsv
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $RNA_MY_ALIGN_DIR/SRR3475377.bam $RNA_MY_REF_GTF > SRR3475377.tsv
```
Subsequently, the tsv files produced by htseq-count are merged with following code and we can see the merged tsv file content
```
join SRR3474779.tsv SRR3474795.tsv | join - SRR3474918.tsv | join - SRR3474972.tsv | join - SRR3475081.tsv | join - SRR3475089.tsv | join - SRR3475326.tsv | join - SRR3475332.tsv | join - SRR3475345.tsv | join - SRR3475362.tsv | join - SRR3475375.tsv | join - SRR3475377.tsv > gene_read_counts_table_all.tsv
echo "GeneID SRR3474779 SRR3474795 SRR3474918 SRR3474972 SRR3475081 SRR3475089 SRR3475326 SRR3475332 SRR3475345 SRR3475362 SRR3475375 SRR3475377" > header.txt
cat header.txt gene_read_counts_table_all.tsv | grep -v "__" | perl -ne 'chomp $_; $_ =~ s/\s+/\t/g; print "$_\n"' > gene_read_counts_table_all_final.tsv
rm -f gene_read_counts_table_all.tsv header.txt
head gene_read_counts_table_all_final.tsv
```

## Differential Expression
Ballgown is going to use for comparison visualizations of the LUSC and normal conditions. First of all, in our R code which is provided in section [TutorialPart1ballgown.R](R_Part1Ballgown.R), we’ll eliminate the low differentially expressed genes with rowVars(bg) >1 function. We use statistical testing to decide whether, for a given gene, an observed difference in read counts is significant, that is, whether it is greater than what would be expected just due to natural random variation. Therefore, in the R codes in section [TutorialPart1ballgown.R](R_Part1Ballgown.R) we only take the genes which are significant(p-value <0.05). Also, we try to eliminate some of the genes by assigning a threshold to FoldChange(FC). FC describes the ratio of two values (ratio of expression in healthy/diseased cases). Due to FC commonly used with log2 if the log2(FC) is equal to 1 between A and B, then A is twice as big as B (or A is 200% of B). Then we’re going to eliminate the genes whether their log2(FC) is greater than 1 or not. If their log2(Foldchange) is less than 1 these genes will be eliminated. But in that case, there would be a biased problem to use FC. For example, if the difference between A and B is 10.000 but the ratio is less than 2 our filtering mechanism will eliminate that gene. On the other hand, if the difference between A and B is 10 but the ratio is 6 our filtering system will account for that gene. Therefore after eliminate the genes by their FC values we’ll sort the genes according to their p-values by increased order and the top 20 genes(lowest p-values) will be printed and the top gene is going to use for Ballgown visualization.

### Differential Expression with HTSEQ-Count
In this section, bioconducter package edgeR have been used for differential expression of count-based RNA-seq data which is generated in section HTSEQ-Count. After running the R codes in [edgeR.R](edgeR.R) section two different significant differentially expressed genes list generated. One of them is generated using Ballgown and the other is generated using edgeR. Below figure provide a nice visualization of overlaps with a venn diagram between Ballgown DE and edgeR DE.

![Ballgown-DE vs edgeR-DE](../main/images/venndiagram.png)

### Visualization of DE Genes with Ballgown
In this section, the summary statistics for FPKM values of genes are plotted and analyzed with Ballgown package and R codes can be provided in [TutorialPart2-ballgown.R](R_Part2Ballgown.R). Also the plots can be provided by follwoing url:[](http://bioinfo05.mu.edu.tr/cihan/Tutorial_Part2_ballgown_output.pdf).

### Analysis and Visualization of DE Genes without Ballgown
In this section, the analysis and visualizations have been progressed without using Ballgown package. All the necessary R codes provided in section [Supplementary.R](Supplementary.R). In this section, 11 different plots are plotted to see discover differences between/among LUSC and Normal samples, examine the differential expression estimates, visualize the expression estimates and highlight those genes that appear to be differentially expressed. All 11 plots can be reached via [](http://bioinfo05.mu.edu.tr/cihan/enrichment_csv_files/Tutorial_Part3_Supplementary_R_output.pdf).

## RNA Sequence Mutation(ANNOVAR)
In this section, de-novo(new) mutation discovery from a given BAM file have been performed. During this section [variation pipeline](http://eng1.mu.edu.tr/~tugba/SeqAnalysis/variation.pipeline) will be followed. Gene-based annotation of genetic variants tries to identify on the human genome, hg38. Gene-based annotation will highlight the exact amino acid change if the mutation is in the exonic region and the predicted effect on the function of the known gene.

## Gene Set Enrichment Analysis(GSEA)
The Gene Set Enrichment Analysis(GSEA) method derives its power by focusing on gene sets, that is, groups of genes that share common biological functions, chromosomal locations, or regulations. In this section three different pathways is going to use for GSEA:
  • Gene Ontology(GO) Pathway
  • KEGG Pathway
  • Immune Pathway
R Codes for Enrichment Analysis provided in section [Gene_Set_Enrichment_Analysis_R](GeneSetEnrinchmentAnalysis.R) Code. Also, necessary output csv files and plot can be observed via [](http://bioinfo05.mu.edu.tr/cihan/GSEA/).






