

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
