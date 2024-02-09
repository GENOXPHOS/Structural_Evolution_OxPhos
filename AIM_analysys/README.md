# Allelic Imbalance (AIM) analysis 
It is required to know the genotypes of analyzed samples, to focus the analysis to variants for which individuals are heterozygotes placed in non overlapping regions between genes. The AIM assessment starts with from bam files coming from scRNAseq data. The analysis is carried out in two steps:

1- Count the reads that present one allele or the other.
2- Contrast of hypothesis of biallelic expression per gene and per cell.
IMPORTANT: In this analysis we assume that we are working with heterozygotes for which we know their genotypes and since the AIM evaluation is at the single cell level, we must be sure that each bam contains single cell information.
## Count reads for each allele.
The reads of each allele are counted using the AIM.read.finder.py script. This script takes as input the directory containing the aligned scRNAseq reads (bam files) for each cell and a VCF file containing the SNPs in heterozygosity that will be considered. The basic command-line for this script is:
```console
$ python AIM.read.finder.py --vcf your_input.vcf --input /your_directory_containing_bam_files/
```
This script will create a new directory called "Output" in your working directory containing txt files with the reads collected from the bam file, matching the SNP positions, annotated with their respective genotypes.
You can check the options of this script by using the --help option:
```console
$python AIM.read.finder.py --help
usage: AIM.read.finder.py [-h] --vcf VCF --input INPUT

options:
  -h, --help            show this help message and exit
  --vcf VCF, -v VCF     Genotype VCF for sample (gzipped with index).
  --input INPUT, -i INPUT
                        Directory containing BAM file with mapped RNA-seq
                        reads
```
## Contrast of hypothesis of biallelic expression per gene and per cell.

This ad hoc method is based on the null hypothesis $H_0$ that genes are biallelic. Then under $H_0$, the probability of expression of each allele is p=0.5. The alternative hypothesis posits that the alleles exhibit unequal probabilities of expression. Specifically, when the gene is monoallelic, resulting in one of the alleles having a significantly higher expression probability than 0.5:

$$
\displaystyle
P(X \leq k | H_0)=
 \sum_{k=1}^n \binom{n}{k} · 0.5 ^ {k} · (1 − 0.5 ) ^ {n − k} 
$$

This analysis is performed by the R script Rscript_4_AIM.R. This Script annotates reads with corresponding mapped genes, perform the AIM analysis per cell and per gene (correcting false discovery rate by Benjamini-Hochberg method, usin a threshold of p<0.05) and create a folder called AIM_output, that contains heatmap representation of AIM analysis and corresponding tsv file with per gene and per cell assessment. In addition, this scrip will evaluate gene-wise trend in all analyzed cells if at least the 75% of cells have a result, terting 3 possible null hypothesis:
1. $H_0$: The count of cells with C57 AIM is equivalent to the count of cells without C57 AIM.
2. $H_0$: The count of cells with Cast AIM is equivalent to the count of cells without Cast AIM.
3. $H_0$: The count of biallelic cells is equivalent to the count of cells without biallelic cells.

In each scenario, when considering the null hypothesis (H0), we assume that the proportion of cells in the analyzed category is equal to p=0.5. However, if the observed proportion significantly exceeds 0.5, with a significance level of 0.05 (p ≤ 0.05), we can infer an alternative hypothesis indicating enrichment within that specific category in our dataset. This script will also output tsv file with gene.wise trend results.
