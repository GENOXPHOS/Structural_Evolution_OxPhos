#######################################################################
## Calculate SEPV as population variability measurement (SEPV.human) ##
#######################################################################


# Load the table for OxPhos subunits annotation
OxPhos=read.delim("/path/to/input_directory/OxPhos.annotation.tsv",sep = "\t",stringsAsFactors = F)

# Define data frame that will contain results from the analysis
Prot.pos = unlist(sapply(1:nrow(OxPhos), function(x){1:OxPhos$length[x]}))
Subunit = unlist(sapply(1:nrow(OxPhos), function(x){rep(OxPhos$Subunit[x],OxPhos$length[x])}))
Gene_Symbol = unlist(sapply(1:nrow(OxPhos), function(x){rep(OxPhos$hgnc_symbol[x],OxPhos$length[x])}))
RC = unlist(sapply(1:nrow(OxPhos), function(x){rep(OxPhos$RC[x],OxPhos$length[x])}))
df=data.frame(Subunit,Gene_Symbol,Prot.pos,RC)



# Define functions to be used:
# Function to determine Shannon Entropy for protein position
ShannonEntropy =function(fqs){-sum(fqs * log2(fqs))}

# Function to determine Shannon Entropy for protein position given an MSA matrix
SEPV_matrix=function(m){return(sapply(1:ncol(m),function(x){ShannonEntropy(table(m[,x])/sum(table(m[,x])))}))}

# Define a function to perform min-max normalization
min_max_norm = function(x) {(x - min(x)) / (max(x) - min(x))}


# Read allele frequencies for missense variants from nuclear OxPhos genes from gnomAD v3.1.2 data 
# These are gene-wise results obtained from gnomAD v.3.1.2 browser
OxPhos.nuclear=read.delim("/path/to/input_directory/gnomAD_v3.1.2_nuclear_OxPhos_subunits.tsv",sep = "\t",stringsAsFactors = F)

# Read allele frequencies for missense variants from mitochondrial OxPhos genes from gnomAD v3.1.2 data 
OxPhos.mt=read.delim("/path/to/input_directory/gnomAD_v3.1.2_mitochondrial_OxPhos_subunits.tsv",sep = "\t",stringsAsFactors = F)
OxPhos.mt=OxPhos.mt[OxPhos.mt$Homoplasmic.Allele.Frequency!=0,]

# Calculate Shannon entropy to define population variability. As nuclear & mt data have different format,SEPV is calculated separately
df$SEPV.human=sapply(paste0(df$Subunit,":",df$Prot.pos),function(x){
  df.sub=OxPhos.nuclear[paste0(OxPhos.nuclear$Ensembl_ID,":",OxPhos.nuclear$Prot.pos)==x,]
  return(ifelse(nrow(df.sub)>0, yes = ShannonEntropy(c(1-sum(df.sub$Allele.Frequency),df.sub$Allele.Frequency)), no=0))
})

# Define mitochondrial subunits as target
mt=df$Subunit %in% OxPhos.mt$Ensembl_ID
df$SEPV.human[mt]=sapply(paste0(df$Subunit[mt],":",df$Prot.pos[mt]),function(x){
  df.sub=OxPhos.mt[paste0(OxPhos.mt$Ensembl_ID,":",OxPhos.mt$Prot.pos)==x,]
  df.sub$Allele.Frequency=df.sub$Homoplasmic.Allele.Frequency
  return(ifelse(nrow(df.sub)>0, yes = ShannonEntropy(c(1-sum(df.sub$Allele.Frequency),df.sub$Allele.Frequency)), no=0))
})

# Min-max normalize SEPV values in human populations
df$SEPV.human=min_max_norm(df$SEPV.human)

# Write table with all results
write.table(df,"/path/to/output_directory/SEPV.human.tsv",sep = "\t",row.names = F,quote = F)


