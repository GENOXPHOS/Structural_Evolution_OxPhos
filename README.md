

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15207045.svg)](https://doi.org/10.5281/zenodo.15207045)


# Structural Diversity and Evolutionary Constrains of Oxidative Phosphorylation

This github repository contains the code used for analysis in the paper entittled"Structural Diversity and Evolutionary Constrains of Oxidative Phosphorylation". Code is provided either in the form of R scripts or Python scripts. All data used for these analyses were included in a zenodo repository (https://doi.org/10.5281/zenodo.15207045).

## Summary of contents
#### OxPhos_variability_1000GP
This folder provides the code used to determine the extent of individual genetic diversity of the OxPhos system subunits in all the 2,504 individualized genotypes included in the 1,000 Genomes Project:

- `Potential_variability_calculation.R` : This script determines the potential nimber of different verssions of each respiratory complex that an individual is able to assemble.

#### SEPV_calculation
This folder provides the code used to determine the extent of individual genetic diversity of the OxPhos system subunits in all the 2,504 individualized genotypes included in the 1,000 Genomes Project:

- `H.from.gnomAD.R` : This R script identifies human constraints at OxPhos subunit positions through Shannon Entropy Position Variation (SEPV), using allele frequency data from gnomAD v3.2.1.
- `H.from.msa.R` : This script computes a conservation score for each position in a multiple sequence alignment (MSA) based on Shannon Entropy Position Variation (SEPV), which quantifies the variability observed at each alignment column.

#### IFRs_detection
This folder provides the code used to identify binding sites in our insilico models: 

- `IFR.R` : This script detects and annotates residues belonging to binding sites in pdb structures and anotate them with additional information.

#### OxPhos_Evolutionary_Strategies_Analyses
This folder provides the code used to explore evolutionary strategies for each OxPhos complexes: 

- `SEPV_analysis1.R` : This script analyzes the evolutionary pathways that have shaped the OxPhos system by comparing Shannon entropy variations across respiratory complexes and subunit types, categorized based on their coding chromosomes (mtDNA, X, or autosomal).
- `IFR_conservation_analysis_1.R` : This script examines the evolutionary patterns of binding sites between interacting subunits, comparing them to other regions within the respiratory complexes (Non-IFRs). The analysis explores whether IFRs have undergone unique evolutionary dynamics compared to the remainder of the protein, while also accounting for the possible impact of subunit type.
- `IFR_conservation_analysis_2.R` : This script captures evolutionary trends from vertebrates to humans through ConScore and examines whether conservation at binding sites differs based on the type of subunits involved, taking into account their unique inheritance mechanisms and evolutionary patterns. The analysis includes a comparison of ConScore behavior against AlphaMissense hotspots.

#### ConScore_benchmarking
This folder provides the code used to evaluate the effectiveness of ConScore compared to AlphaMissense hotspots, this analysis focuses on missense variants in ClinVar impacting the OxPhos system. Performance is assessed using the area under the receiver operating characteristic curve (ROC-AUC) and the area under the precision-recall curve (PR-AUC): 

- `ConScore_benchmarking.R` : This script is used for benchmarking ConScore.

#### PSI_analysis
This folder provides the code and files used to analyze the relationship between population variability, level of conservation and impact on structure:

- `PSI_calculation.R` : This script assesses the resilience of proteins to mutations by evaluating the stability of the protein structure at each subunit position.
- `PSI.analysis.R` : This script examines the interplay between population variability, conservation levels, and structural impact.

#### Extreme_Genetic_Constriction_Strategies
This folder provides the code and files used to study extreme strategies for shaping human population constraints in OxPhos system:

- `Extreme_Strategies_analysis_1.R` : This script compiles positions exhibiting extreme constraint strategies (homozygous-biased and fixed positions) in human populations.
- `Extreme_Strategies_analysis_2.R` : This script evaluates if fixed positions are randomly distributed among RCs or subunits.

#### AIM_analysis
This folder provides the code used to explore Allelic Imbalance (AIM) as a phenomenon capable of influencing the effective OxPhos variability that is finally assembled:

- `AIM.read.finder.py` : This script counts the reads that present one allele or the other.
- `Rscript_4_AIM.R` : This script performs contrast of hypothesis of biallelic expression per gene and per cell.


## How to cite

If you make use of analysis scripts or data from this work, please cite as follows:

Not available

## Contact: 
jaenriquez@cnic.es or jlcabreraa@cnic.es 

## References:

- Knaus BJ, Grünwald NJ. vcfr: a package to manipulate and visualize variant call format data in R. Mol Ecol Resour. 2017 Jan;17(1):44-53. doi: 10.1111/1755-0998.12549. Epub 2016 Jul 12. PMID: 27401132.

- Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York. ISBN 978-3-319-24277-4, https://ggplot2.tidyverse.org.

- Wickham H, François R, Henry L, Müller K, Vaughan D (2023). dplyr: A Grammar of Data Manipulation. R package version 1.1.4, https://github.com/tidyverse/dplyr, https://dplyr.tidyverse.org.

- Charif D, Lobry J (2007). “SeqinR 1.0-2: a contributed package to the R project for statistical computing devoted to biological sequences retrieval and analysis.” In Bastolla U, Porto M, Roman H, Vendruscolo M (eds.), Structural approaches to sequence evolution: Molecules, networks, populations, series Biological and Medical Physics, Biomedical Engineering, 207-232. Springer Verlag, New York. ISBN : 978-3-540-35305-8.

- Pagès H, Aboyoun P, Gentleman R, DebRoy S (2024). Biostrings: Efficient manipulation of biological strings. doi:10.18129/B9.bioc.Biostrings <https://doi.org/10.18129/B9.bioc.Biostrings>, R package version 2.72.1, <https://bioconductor.org/packages/Biostrings>.

- Kassambara A (2023). ggpubr: 'ggplot2' Based Publication Ready Plots. R package version 0.6.0, https://rpkgs.datanovia.com/ggpubr/.

- Ogle DH, Doll JC, Wheeler AP, Dinno A (2025). FSA: Simple Fisheries Stock Assessment Methods. R package version 0.9.6, https://CRAN.R-project.org/package=FSA.

- Grant BJ, Rodrigues AP, ElSawy KM, McCammon JA, Caves LS. Bio3d: an R package for the comparative analysis of protein structures. Bioinformatics. 2006 Nov 1;22(21):2695-6. doi: 10.1093/bioinformatics/btl461. Epub 2006 Aug 29. PMID: 16940322.

- R Core Team (2024). R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria. <https://www.R-project.org/>.

- Saito T, Rehmsmeier M. Precrec: fast and accurate precision-recall and ROC curve calculations in R. Bioinformatics. 2017 Jan 1;33(1):145-147. doi: 10.1093/bioinformatics/btw570. Epub 2016 Sep 1. PMID: 27591081; PMCID: PMC5408773.

- Sing T, Sander O, Beerenwinkel N, Lengauer T. ROCR: visualizing classifier performance in R. Bioinformatics. 2005 Oct 15;21(20):3940-1. doi: 10.1093/bioinformatics/bti623. Epub 2005 Aug 11. PMID: 16096348.

- Ben-Shachar et al., (2020). effectsize: Estimation of Effect Size Indices and Standardized Parameters. Journal of Open Source Software, 5(56), 2815, https://doi.org/10.21105/joss.02815.

- Python Software Foundation. (2023). os — Miscellaneous operating system interfaces. En Python Documentation. Recuperado de https://docs.python.org/3/library/os.html

- Python Software Foundation. (2023). argparse — Parser for command-line options, arguments and sub-commands. En Python Documentation. Recuperado de https://docs.python.org/3/library/argparse.html

- Heger A, Pysam developers. pysam: a Python module for reading and manipulating SAM/BAM/VCF/BCF files. Available at: https://pysam.readthedocs.io.

- McKinney, W. (2010). Data Structures for Statistical Computing in Python. Proceedings of the 9th Python in Science Conference, 51-56.

- Gu Z, Eils R, Schlesner M. Complex heatmaps reveal patterns and correlations in multidimensional genomic data. Bioinformatics. 2016 Sep 15;32(18):2847-9. doi: 10.1093/bioinformatics/btw313. Epub 2016 May 20. PMID: 27207943.

- Gu Z, Gu L, Eils R, Schlesner M, Brors B. circlize Implements and enhances circular visualization in R. Bioinformatics. 2014 Oct;30(19):2811-2. doi: 10.1093/bioinformatics/btu393. Epub 2014 Jun 14. PMID: 24930139.

- Davis T (2024). optparse: Command Line Option Parser. R package version 1.7.5, <https://CRAN.R-project.org/package=optparse>.

- Barrett T, Dowle M, Srinivasan A, Gorecki J, Chirico M, Hocking T, Schwendinger B (2024). data.table: Extension of `data.frame`. R package version 1.16.2, <https://CRAN.R-project.org/package=data.table>.


