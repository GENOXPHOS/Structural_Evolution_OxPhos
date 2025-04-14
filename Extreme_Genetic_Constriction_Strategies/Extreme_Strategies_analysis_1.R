library("dplyr")

#######################################################################
## Calculate SEPV as population variability measurement (SEPV.human) ##
#######################################################################


# Load the table for OxPhos subunits annotation
# OxPhos=read.delim("/path/to/input_directory/OxPhos.annotation.tsv",sep = "\t",stringsAsFactors = F)
OxPhos=read.delim("/path/to/input_directory/OxPhos.annotation.tsv",sep = "\t",stringsAsFactors = F)

# Define data frame that will contain results from the analysis
Prot.pos = unlist(sapply(1:nrow(OxPhos), function(x){1:OxPhos$length[x]}))
Subunit = unlist(sapply(1:nrow(OxPhos), function(x){rep(OxPhos$ensembl_gene_id[x],OxPhos$length[x])}))
Gene_Symbol = unlist(sapply(1:nrow(OxPhos), function(x){rep(OxPhos$hgnc_symbol[x],OxPhos$length[x])}))
RC = unlist(sapply(1:nrow(OxPhos), function(x){rep(OxPhos$RC[x],OxPhos$length[x])}))
df=data.frame(Subunit,Gene_Symbol,Prot.pos,RC)


# Define functions to be used

# Function to test the number of observed Homozygotes is equal to the number of expected homozygotes
HBP=function(df.input){
  N=round(x = max(df.input$Allele.Number)/2,digits = 0)
  Obs.alt.homo=sum(df.input$Homozygote.Count)
  Obs.alt.het=sum(df.input$Allele.Count-df.input$Homozygote.Count*2)
  Obs.ref.homo = N-(Obs.alt.homo+Obs.alt.het)
  Ner_Homo.obs = Obs.ref.homo + Obs.alt.homo
  Ref.allele.fq = (1-sum(df.input$Allele.Frequency))
  Alt.allele.fq = df.input$Allele.Frequency
  Expected_homo=sum(c(Ref.allele.fq,Alt.allele.fq)**2)
  return(binom.test(x = Ner_Homo.obs, n = N, p=Expected_homo,alternative = 'greater')$p.value)
}

# Function to determine Shannon Entropy for protein position
ShannonEntropy =function(fqs){-sum(fqs * log2(fqs))}

# Function to determine Shannon Entropy for protein position given an MSA matrix
SEPV_matrix=function(m){return(sapply(1:ncol(m),function(x){ShannonEntropy(table(m[,x])/sum(table(m[,x])))}))}

# Define a function to perform min-max normalization
min_max_norm = function(x) {(x - min(x)) / (max(x) - min(x))}

# Read allele frequencies for missense variants from nuclear OxPhos genes from gnomAD v3.1.2 data 
# These are gene-wise results obtained from gnomAD v.3.1.2 browser
# OxPhos.nuclear=read.delim("/path/to/input_directory/gnomAD_v3.1.2_nuclear_OxPhos_subunits.tsv",sep = "\t",stringsAsFactors = F)
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


##############################################
## Identification of Non changing positions ##
##############################################

# Define in which protein positions there are not variants described in gnomAD v3.1.2
df$Fixed=ifelse(df$SEPV.human>0, yes = 0, no=1)


###################################################
## Identification of Homozygous biased positions ##
###################################################


# Define mitochondrial and X-linked OxPhox genes 
Mt.genes=c("ENSG00000198888", "ENSG00000198763", "ENSG00000198804", "ENSG00000198712", "ENSG00000228253", "ENSG00000198899", "ENSG00000198938", "ENSG00000198840", "ENSG00000212907","ENSG00000198886", "ENSG00000198786", "ENSG00000198695", "ENSG00000198727")
X.genes=c("ENSG00000125356", "ENSG00000131174", "ENSG00000147123")

# Set filter for autosomal genes to analyze the existence of homozygous biased positions
f1 = df$Subunit %in% X.genes == F
f2 = paste0(df$Subunit,':',df$Prot.pos) %in% paste0(OxPhos.nuclear$Ensembl_ID,':',OxPhos.nuclear$Prot.pos)

# Assess hozygous biased positions to limit population variability as binomial test
df$p.value[f1 & f2]=sapply(paste0(df$Subunit[f1 & f2],':',df$Prot.pos[f1 & f2]),function(x){
  df.sub=OxPhos.nuclear[paste0(OxPhos.nuclear$Ensembl_ID,':',OxPhos.nuclear$Prot.pos)==x,]
  return(ifelse(nrow(df.sub)>0, yes = HBP(df.sub), no=NA))
})

# Correct False discovery rate by Benjamini-Hochberg method
df$adj.p=p.adjust(df$p.value, method = "BH")

# Classify protein positions with heterozygous penalyzation
df$`Homozygous biased Position`=0
df$`Homozygous biased Position`[df$adj.p<=0.05]=1

df=df %>% group_by(Subunit) %>% summarise(Gene_Symbol=first(Gene_Symbol), RC=first(RC),Fixed=sum(Fixed), `Homozygous biased Position` = sum(`Homozygous biased Position`))

# Write table with all results
write.table(df,"/path/to/output_directory/Table_S5.tsv",sep = "\t",row.names = F,quote = F)


