library(dplyr)
library(seqinr)
library(Biostrings)

#############################################
### Calculate SEPV from MUSCLE alignments ###
#############################################

# Define the directory path where the MUSCLE output files are located
path.muscle="/path/to/input_directory/"

# List all files in the given directory that have the ".out" extension.
# NOTE: ".out" is the custom extension used for MUSCLE alignment output files,
# and each file is named using its corresponding Ensembl ID.

l=list.files(path = path.muscle,pattern = ".out")

# Create an empty data frame to store SEPV data for all OxPhos subunits
df=data.frame()

# Define a function to perform min-max normalization
min_max_norm = function(x) {(x - min(x)) / (max(x) - min(x))}

# Define a function to calculate Shannon entropy from a vector of frequencies

ShannonEntropy =function(vec) {-sum(vec * log2(vec))}

# Define a function called "wtable" that computes weighted (adjusted) frequencies 
# used in the calculation of Shannon entropy.

# The function accepts two parameters: "x", representing the observed amino acid 
# frequencies at a specific position, and "w", a set of weights applied to adjust 
# these frequencies based on how dissimilar each sequence is from the reference 
# (typically the human sequence, which is the first in the alignment).

# This weighting correction ensures that entropy values are comparable across 
# different proteins, even when the corresponding sequence for a particular species 
# is missing in some alignments.

wtable=function(x,w){
  v1=x
  v2=sapply(unique(v1), function(aa){return(sum(w[which(v1==aa)]) / sum(w))})
  return(v2)
}

# Load table for dataframe needed for subunits Annotation 
OxPhos=read.delim("/path/to/input_directory/OxPhos.annotation.tsv",sep = "\t",stringsAsFactors = F)

for (f in l){
  
  
  # Define Subunit name as Ensembl ID
  Gene.ID=gsub(".out","",fixed = T,f)
  
  # Load MSA 
  msa = read.alignment(file = paste0(path.muscle,f), format = "fasta")
  
  # Calculate weights for frequencies correction based on how dissimilar each sequence is from 
  # the reference (allways the first sequence)
  seqs = AAStringSet(sapply(msa$seq, as.character))
  dist_matrix = as.matrix(dist.alignment(msa, matrix = "identity",gap=T))
  w=1-dist_matrix[,1]
  msa.matrix=as.matrix.alignment(msa)
  
  # We take as reference the human sequence which is the first one to fix the protein positions
  msa.matrix=msa.matrix[,msa.matrix[1,]!="-"]
  colnames(msa.matrix)=c(1:ncol(msa.matrix))
  
  # Calculate SEPV by protein position
  Prot.pos=c()
  SEPV=c()
  for (c in 1:ncol(msa.matrix)){
    Prot.pos=c(Prot.pos,c)
    SEPV=c(SEPV,ShannonEntropy(wtable(as.vector(msa.matrix[,c]),w)))
  }
  
  # Create a data frame "df1" with SEPV information for working subunit with 
  df1=data.frame(Prot.pos=Prot.pos,SEPV=SEPV)
  
  # Add identification details to the subunit dataframe
  df1$Subunit=Gene.ID
  df1$entrezgene_id=OxPhos$entrezgene_id[OxPhos$ensembl_gene_id==Gene.ID]
  df1$Gene_symbol=OxPhos$hgnc_symbol[OxPhos$ensembl_gene_id==Gene.ID]
  
  # Order df variables
  df1=df1[,c("Subunit","Gene_symbol","entrezgene_id","Prot.pos","SEPV")]
  
  # Update output with subunit information 
  df=rbind(df,df1)
}

# Min-max normalize SEPV values
df$SEPV=min_max_norm(df$SEPV)

# Write output file in tsv format
write.table(df,"/path/to/output_directory/SEPV.output.tsv",sep = "\t",row.names = F,quote = F)
