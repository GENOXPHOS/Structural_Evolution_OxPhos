library(dplyr)

# Gather energy variation data for each subunit analyzed using a PyRosetta script (predict_ompLA_ddG.py),
# customized for our purposes, into a dataframe object from the DDG directory
path="/path/to/input_directory/DDG_data"
l=list.files(path = path,pattern = "_ddG.out")

df=data.frame()
for (f in l){
  df.1=read.delim(paste0(path,"/",f),sep = " ",stringsAsFactors = F)
  gene=gsub("_ddG.out","",f)
  df.1=df.1[which(df.1$Ref!=df.1$AA),]
  df.1=df.1[which(df.1$Ref!="Z"),]
  df.1$AA.sub=paste0(df.1$Ref,df.1$AA)
  df.1$Gene=gene
  df=rbind(df,df.1)
}
df$Gene=gsub("MT_","MT-",df$Gene)


# Function to calculate the p-value for the Kolmogorov-Smirnov test, representing the probability 
# that the impact on protein structure stability, caused by the mutation of an AA at a specific 
# position to one of the 19 alternative AAs, exceeds the impact observed from mutating all protein AAs 
# to the remaining 19 alternatives
calculate_ks_pvalue = function(ddG_group, ddG_all) {
  suppressWarnings(ks.test(abs(ddG_group), abs(ddG_all), alternative = "less", exact = FALSE)$p.value)
}

# Group by Gene and Res, then compute the p-value for each group
df = df %>%
  group_by(Gene, Res) %>%
  mutate(PSI = calculate_ks_pvalue(df$ddG, ddG)) %>%
  ungroup()
end.time = Sys.time()
time.taken = end.time - start.time
time.taken

# Collapse the dataframe df to retain a single value per position
df=unique(df[,c("Gene","Res","PSI")])

# Load table for dataframe needed for subunits Annotation 
OxPhos=read.delim("./OxPhox.tsv", stringsAsFactors = F)

# Annotate dataframe 
for (g in unique(df$Gene)){
  df$Ensembl_ID[df$Gene==g]=rep(OxPhos$ensembl_gene_id[OxPhos$hgnc_symbol==g], length(which(df$Gene==g)))
  df$RC[df$Gene==g]=rep(OxPhos$RC[OxPhos$hgnc_symbol==g], length(which(df$Gene==g)))
}

# Save results of Probability of Structural Impact (PSI) as tsv file, reordering columns
write.table(df[,c("Ensembl_ID","Gene","Res","PSI","RC")],"/path/to/output_directory/Table_S4.tsv",sep = "\t",quote = F,row.names = F)
