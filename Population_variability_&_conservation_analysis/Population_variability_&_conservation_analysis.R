# Load libraries
library(parallel)
library(ggplot2)
library(dplyr)
library(seqinr)
library(bio3d)

##################################################################################################################
## This script performs all the analyses considering that we compare values that are linked to residues that    ##
## will compose different multiprotein complexes, where some of the subunits that compose them will be repeated ##
## or have alternative isoforms encoded by paralogs. Consequently, when comparing the different distributions   ##
## in each hypothesis test, the stoichiometric composition of the respiratory complexes was taken into account. ##                                                                                     ##
##################################################################################################################


# Load table for OxPhos Annotation. This table used to have a unified annotation with the gene symbol, specially for CV subunits
OxPhos=read.delim('./OxPhos.annotation.tsv',sep = '\t',stringsAsFactors = F)

# Define data frame that will contain results from the analysis
Prot.pos = unlist(sapply(1:nrow(OxPhos), function(x){1:OxPhos$length[x]}))
Subunit = unlist(sapply(1:nrow(OxPhos), function(x){rep(OxPhos$Subunit[x],OxPhos$length[x])}))
Gene_Symbol = unlist(sapply(1:nrow(OxPhos), function(x){rep(OxPhos$hgnc_symbol[x],OxPhos$length[x])}))
RC = unlist(sapply(1:nrow(OxPhos), function(x){rep(OxPhos$RC[x],OxPhos$length[x])}))
df=data.frame(Subunit,Gene_Symbol,Prot.pos,RC)

# Create output directories
dir.create('./Output')
dir.create('./Output/IFRs_SEPV_analysis')
dir.create('./Output/Whole_RC_SEPV_analysis')

# Define mitochondrial and X-linked OxPhox genes 
Mt.genes=c("ENSG00000198888", "ENSG00000198763", "ENSG00000198804", "ENSG00000198712", "ENSG00000228253", "ENSG00000198899", "ENSG00000198938", "ENSG00000198840", "ENSG00000212907","ENSG00000198886", "ENSG00000198786", "ENSG00000198695", "ENSG00000198727")
X.genes=c("ENSG00000125356", "ENSG00000131174", "ENSG00000147123")

################################################################################
## Calculate SEPV as population variability measurement (SEPV.gnomAD v-3.1.2) ##
################################################################################

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

# Read allele frequencies for missense variants from nuclear OxPhos genes from gnomAD v3.1.2 data 
# These are gene-wise results obtained from gnomAD v.3.1.2 browser
OxPhos.nuclear=read.delim('./gnomAD_v3.1.2_nuclear_OxPhos_subunits.tsv',sep = '\t',stringsAsFactors = F)

# Read allele frequencies for missense variants from mitochondrial OxPhos genes from gnomAD v3.1.2 data 
OxPhos.mt=read.delim('./gnomAD_v3.1.2_mitochondrial_OxPhos_subunits.tsv',sep = '\t',stringsAsFactors = F)

# Calculate Shannon entropy to define population variability. As nuclear & mt data have different format,SEPV is calculated separately
df$SEPV.gnomAD.v3=sapply(paste0(df$Subunit,':',df$Prot.pos),function(x){
  df.sub=OxPhos.nuclear[paste0(OxPhos.nuclear$Ensembl_ID,':',OxPhos.nuclear$Prot.pos)==x,]
  return(ifelse(nrow(df.sub)>0, yes = ShannonEntropy(c(1-sum(df.sub$Allele.Frequency),df.sub$Allele.Frequency)), no=0))
})

# Define mitochondrial subunits as target
mt=df$Subunit %in% OxPhos.mt$Ensembl_ID
df$SEPV.gnomAD.v3[mt]=sapply(paste0(df$Subunit[mt],':',df$Prot.pos[mt]),function(x){
  df.sub=OxPhos.mt[paste0(OxPhos.mt$Ensembl_ID,':',OxPhos.mt$Prot.pos)==x,]
  df.sub$Allele.Frequency=(df.sub$Homoplasmic.Allele.Count+df.sub$Heteroplasmic.Allele.Count)/df.sub$Allele.Number
  return(ifelse(nrow(df.sub)>0, yes = ShannonEntropy(c(1-sum(df.sub$Allele.Frequency),df.sub$Allele.Frequency)), no=0))
})


##############################################
## Identification of Non changing positions ##
##############################################

# Define in which protein positions there are not variants described in gnomAD v3.1.2
df$Non.Changing.Position=ifelse(df$SEPV.gnomAD.v3>0, yes = 0, no=1)

###################################################
## Identification of Homozygous biased positions ##
###################################################

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
df$Homozygous.biased.Position=0
df$Homozygous.biased.Position[df$adj.p<=0.05]=1


##################################################################################
## Calculate SEPV as conservation measurement from MSA (SEPV.msa.interespecies) ##
##################################################################################

l=list.files(path = './muscle/',pattern = '.out')
l=l[which(gsub('.out','',l) %in% unique(df$Subunit))]

# Calculate SEPV from MSA taking as reference the first species that in all the files must be the human sequence
df$SEPV.msa.interspecies=unlist(sapply(l,function(x){
  gene=gsub('.out','',x)
  msa <- read.alignment(file = paste0('./muscle/',x), format = "fasta")
  msa.matrix=as.matrix.alignment(msa)
  msa.matrix=msa.matrix[,msa.matrix[1,]!='-']
  return(SEPV_matrix(msa.matrix))
}))

# Write table with all results
write.table(df,'./Output/Table_S2.tsv',sep = '\t',row.names = F,quote = F)

# Summarize of NCP and HBP results
df2=df[,c("Subunit","Gene_Symbol","RC","Non.Changing.Position","Homozygous.biased.Position")] %>% group_by(Subunit,Gene_Symbol,RC) %>% summarise(NCP = sum(Non.Changing.Position, na.rm = TRUE),HBP=sum(Homozygous.biased.Position, na.rm = TRUE),length=n())

# Estimate the number of non changing positions by protein size
df2$normalized.NCP=df2$NCP/df2$length

# Write table with gene-wise summary of NCP and HBP results
write.table(df2,'./Output/Table_S1.tsv',sep = '\t',row.names = F,quote = F)


# Set proper stoichiometry for RCs regarding their composition of subunits
df2=rbind(df2,df2[df2$Subunit=='ENSG00000004779',], df2[rep(rownames(df2[df2$Subunit=='ENSG00000152234',]),2),], df2[df2$RC=='CIII',],
         df2[rep(rownames(df2[df2$Subunit=='ENSG00000110955',]),2),], df2[rep(rownames(df2[df2$Subunit=='ENSG00000159199',]),7),], 
         df2[rep(rownames(df2[df2$Subunit=='ENSG00000135390',]),7),], df2[rep(rownames(df2[df2$Subunit=='ENSG00000154518',]),7),])


# Define a function to obtain a data frame with 3 columns: Conditions compared, Alternative hypothesis and Wilcox test p.value
wilcox_test_func <- function(vector1, vector2) {
  
  # Perform Wilcoxon test for 'less'
  test_less <- wilcox.test(vector1, vector2, alternative = "less")
  
  # Perform Wilcoxon test for 'greater'
  test_greater <- wilcox.test(vector1, vector2, alternative = "greater")
  
  # Check if the p-value for 'less' is significant
  if (test_less$p.value < 0.05) {
    results <- data.frame(Result=paste0(deparse(substitute(vector1)),' Vs ',deparse(substitute(vector2))),p_value = test_less$p.value, hypothesis = "less")
  }
  
  # Check if the p-value for 'greater' is significant
  if (test_greater$p.value < 0.05) {
    results <- data.frame(Result=paste0(deparse(substitute(vector1)),' Vs ',deparse(substitute(vector2))), p_value = test_greater$p.value, hypothesis = "greater")
  }
  
  # If no significant hypothesis was found, return 'ns'
  if (test_greater$p.value > 0.05 & test_less$p.value > 0.05) {
    results <- data.frame(Result= paste0(deparse(substitute(vector1)),' Vs ',deparse(substitute(vector2))),p_value = NA, hypothesis = "ns")
  }
  
  return(results)
}

# Set types of subunits by chromosome/inheritance/evolutionary strategy

df$Chromosome="Autosomal"
df$Chromosome[df$Subunit %in% Mt.genes]="Mitochondrial"
df$Chromosome[df$Subunit %in% X.genes]="X"

df2$Genome="Autosomal"
df2$Genome[df2$Subunit %in% Mt.genes]="Mitochondrial"
df2$Genome[df2$Subunit %in% X.genes]="X"


# Define distributions to be compared by RC for the number of Non changing positions normalized by size
NNCP.CI=df2$normalized.NCP[df2$RC=='CI']
NNCP.CII=df2$normalized.NCP[df2$RC=='CII']
NNCP.CIII=df2$normalized.NCP[df2$RC=='CIII']
NNCP.CIV=df2$normalized.NCP[df2$RC=='CIV']
NNCP.CV=df2$normalized.NCP[df2$RC=='CV']

# Create a data frame with all 10 comparing the number of non changing positions normalized by protein size between RCs
results=data.frame()
results=rbind(results,wilcox_test_func(NNCP.CI,NNCP.CII))
results=rbind(results,wilcox_test_func(NNCP.CI,NNCP.CIII))
results=rbind(results,wilcox_test_func(NNCP.CI,NNCP.CIV))
results=rbind(results,wilcox_test_func(NNCP.CI,NNCP.CV))
results=rbind(results,wilcox_test_func(NNCP.CII,NNCP.CIII))
results=rbind(results,wilcox_test_func(NNCP.CII,NNCP.CIV))
results=rbind(results,wilcox_test_func(NNCP.CII,NNCP.CV))
results=rbind(results,wilcox_test_func(NNCP.CIII,NNCP.CIV))
results=rbind(results,wilcox_test_func(NNCP.CIII,NNCP.CV))
results=rbind(results,wilcox_test_func(NNCP.CIV,NNCP.CV))

# Write tsv table with the results
write.table(results,'./Output/Table_S3.tsv',sep = '\t',quote = F,row.names = F)

# Boxplot comparing the number of Non Changing Positions normalized by subunit size 
pdf('./Output/NNCP.by.RC.pdf',width = 8,height = 6)
ggplot(df2,aes(x=RC, y=normalized.NCP,fill=RC))+
  geom_boxplot(alpha=0.8, outlier.colour="red",outlier.fill="red") + 
  ylab("Ner of Non Changing Positions normalized by protein size") + 
  xlab('') + 
  theme(legend.position = "none")
dev.off()


####################################################################
#### Compare SEPV.gnomAD.v3 & SEPV.msa.interspecies between RCs ####
####################################################################

# Define values of both meanings SEPV by RC considering stoichiometric presence of subunits in the structure

# SEPV (both SEPV.gnomAD.v3 & SEPV.msa.interspecies) in CI 
CI.SEPV.gnomAD.v3=c(df$SEPV.gnomAD.v3[which(df$RC=='CI')],df$SEPV.gnomAD.v3[df$Subunit=='ENSG00000004779'])
CI.SEPV.msa.interspecies=c(df$SEPV.msa.interspecies[which(df$RC=='CI')],df$SEPV.msa.interspecies[df$Subunit=='ENSG00000004779'])

# SEPV (both SEPV.gnomAD.v3 & SEPV.msa.interspecies) in CII
CII.SEPV.gnomAD.v3=df$SEPV.gnomAD.v3[which(df$RC=='CII')]
CII.SEPV.msa.interspecies=df$SEPV.msa.interspecies[which(df$RC=='CII')]

# SEPV (both SEPV.gnomAD.v3 & SEPV.msa.interspecies) in CIII
CIII.SEPV.gnomAD.v3=c(df$SEPV.gnomAD.v3[which(df$RC=='CIII')],df$SEPV.gnomAD.v3[which(df$RC=='CIII')])
CIII.SEPV.msa.interspecies=c(df$SEPV.msa.interspecies[which(df$RC=='CIII')],df$SEPV[which(df$RC=='CIII')])

# SEPV (both SEPV.gnomAD.v3 & SEPV.msa.interspecies) in CIV
CIV.SEPV.gnomAD.v3=df$SEPV.gnomAD.v3[which(df$RC=='CIV')]
CIV.SEPV.msa.interspecies=df$SEPV.msa.interspecies[which(df$RC=='CIV')]

# SEPV (both SEPV.gnomAD.v3 & SEPV.msa.interspecies) in CV 
CV.SEPV.gnomAD.v3=c(df$SEPV.gnomAD.v3[which(df$RC=='CV')],rep(df$SEPV.gnomAD.v3[df$Subunit=='ENSG00000152234'],2),
                    rep(df$SEPV.gnomAD.v3[df$Subunit=='ENSG00000110955'],2), rep(df$SEPV.gnomAD.v3[df$Subunit=='ENSG00000159199'],7), 
                    rep(df$SEPV.gnomAD.v3[df$Subunit=='ENSG00000135390'],7), rep(df$SEPV.gnomAD.v3[df$Subunit=='ENSG00000154518'],7))

CV.SEPV.msa.interspecies=c(df$SEPV.msa.interspecies[which(df$RC=='CV')],rep(df$SEPV.msa.interspecies[df$Subunit=='ENSG00000152234'],2),
                           rep(df$SEPV.msa.interspecies[df$Subunit=='ENSG00000110955'],2), rep(df$SEPV.msa.interspecies[df$Subunit=='ENSG00000159199'],7), 
                           rep(df$SEPV.msa.interspecies[df$Subunit=='ENSG00000135390'],7), rep(df$SEPV.msa.interspecies[df$Subunit=='ENSG00000154518'],7))


# Create a data frame with all 10 possible comparisons between RCs for SEPV msa interspecies
results=data.frame()
results=rbind(results,wilcox_test_func(CI.SEPV.msa.interspecies,CII.SEPV.msa.interspecies))
results=rbind(results,wilcox_test_func(CI.SEPV.msa.interspecies,CIII.SEPV.msa.interspecies))
results=rbind(results,wilcox_test_func(CI.SEPV.msa.interspecies,CIV.SEPV.msa.interspecies))
results=rbind(results,wilcox_test_func(CI.SEPV.msa.interspecies,CV.SEPV.msa.interspecies))
results=rbind(results,wilcox_test_func(CII.SEPV.msa.interspecies,CIII.SEPV.msa.interspecies))
results=rbind(results,wilcox_test_func(CII.SEPV.msa.interspecies,CIV.SEPV.msa.interspecies))
results=rbind(results,wilcox_test_func(CII.SEPV.msa.interspecies,CV.SEPV.msa.interspecies))
results=rbind(results,wilcox_test_func(CIII.SEPV.msa.interspecies,CIV.SEPV.msa.interspecies))
results=rbind(results,wilcox_test_func(CIII.SEPV.msa.interspecies,CV.SEPV.msa.interspecies))
results=rbind(results,wilcox_test_func(CIV.SEPV.msa.interspecies,CV.SEPV.msa.interspecies))

# Write tsv table with the results
write.table(results,'./Output/Whole_RC_SEPV_analysis/SEPV.msa.interspecies.between.RCs.tsv',sep = '\t',quote = F,row.names = F)

# Create a data frame with all 10 possible comparisons between RCs for SEPV gnomAD v3 values
results=data.frame()
results=rbind(results,wilcox_test_func(CI.SEPV.gnomAD.v3,CII.SEPV.gnomAD.v3))
results=rbind(results,wilcox_test_func(CI.SEPV.gnomAD.v3,CIII.SEPV.gnomAD.v3))
results=rbind(results,wilcox_test_func(CI.SEPV.gnomAD.v3,CIV.SEPV.gnomAD.v3))
results=rbind(results,wilcox_test_func(CI.SEPV.gnomAD.v3,CV.SEPV.gnomAD.v3))
results=rbind(results,wilcox_test_func(CII.SEPV.gnomAD.v3,CIII.SEPV.gnomAD.v3))
results=rbind(results,wilcox_test_func(CII.SEPV.gnomAD.v3,CIV.SEPV.gnomAD.v3))
results=rbind(results,wilcox_test_func(CII.SEPV.gnomAD.v3,CV.SEPV.gnomAD.v3))
results=rbind(results,wilcox_test_func(CIII.SEPV.gnomAD.v3,CIV.SEPV.gnomAD.v3))
results=rbind(results,wilcox_test_func(CIII.SEPV.gnomAD.v3,CV.SEPV.gnomAD.v3))
results=rbind(results,wilcox_test_func(CIV.SEPV.gnomAD.v3,CV.SEPV.gnomAD.v3))

# Write tsv table with the results
write.table(results,'./Output/Whole_RC_SEPV_analysis/SEPV.gnomAD.v3.between.RCs.tsv',sep = '\t',quote = F,row.names = F)

################################################################
#### Compare SEPV (both meanings) between types of subunits ####
################################################################

# SEPV (both SEPV.gnomAD.v3 & SEPV.msa.interspecies) in autosomal 
SEPV.gnomAD.v3.Autosomal=c(df$SEPV.gnomAD.v3[which(df$Chromosome=='Autosomal')],df$SEPV.gnomAD.v3[df$Subunit=='ENSG00000004779'],rep(df$SEPV.gnomAD.v3[df$Subunit=='ENSG00000152234'],2),
                           rep(df$SEPV.gnomAD.v3[df$Subunit=='ENSG00000110955'],2), rep(df$SEPV.gnomAD.v3[df$Subunit=='ENSG00000159199'],7), 
                           rep(df$SEPV.gnomAD.v3[df$Subunit=='ENSG00000135390'],7), rep(df$SEPV.gnomAD.v3[df$Subunit=='ENSG00000154518'],7),
                           df$SEPV.gnomAD.v3[which(df$Chromosome=='Autosomal' & df$RC=='CIII')])
SEPV.msa.interspecies.Autosomal=c(df$SEPV.msa.interspecies[which(df$Chromosome=='Autosomal')],df$SEPV.msa.interspecies[df$Subunit=='ENSG00000004779'],rep(df$SEPV.msa.interspecies[df$Subunit=='ENSG00000152234'],2),
                                  rep(df$SEPV.msa.interspecies[df$Subunit=='ENSG00000110955'],2), rep(df$SEPV.msa.interspecies[df$Subunit=='ENSG00000159199'],7), 
                                  rep(df$SEPV.msa.interspecies[df$Subunit=='ENSG00000135390'],7), rep(df$SEPV.msa.interspecies[df$Subunit=='ENSG00000154518'],7),
                                  df$SEPV.msa.interspecies[which(df$Chromosome=='Autosomal' & df$RC=='CIII')])

# SEPV (both SEPV.gnomAD.v3 & SEPV.msa.interspecies) in Mitochondrial 
SEPV.gnomAD.v3.Mitochondrial=c(df$SEPV.gnomAD.v3[which(df$Chromosome=='Mitochondrial')],df$SEPV.gnomAD.v3[which(df$Chromosome=='Mitochondrial' & df$RC=='CIII')])
SEPV.msa.interspecies.Mitochondrial=c(df$SEPV.msa.interspecies[which(df$Chromosome=='Mitochondrial')],df$SEPV.msa.interspecies[which(df$Chromosome=='Mitochondrial' & df$RC=='CIII')])

# SEPV (both SEPV.gnomAD.v3 & SEPV.msa.interspecies) in X-linked
SEPV.gnomAD.v3.X=df$SEPV.gnomAD.v3[which(df$Chromosome=='X')]
SEPV.msa.interspecies.X=df$SEPV.msa.interspecies[which(df$Chromosome=='X')]

# Create a data frame with all 3 possible comparisons between types of subunits for SEPV msa interspecies
results=data.frame()
results=rbind(results,wilcox_test_func(SEPV.msa.interspecies.Autosomal,SEPV.msa.interspecies.Mitochondrial))
results=rbind(results,wilcox_test_func(SEPV.msa.interspecies.Autosomal,SEPV.msa.interspecies.X))
results=rbind(results,wilcox_test_func(SEPV.msa.interspecies.Mitochondrial,SEPV.msa.interspecies.X))

# Write tsv table with the results
write.table(results,'./Output/Whole_RC_SEPV_analysis/SEPV.msa.interspecies.between.subunit.type.tsv',sep = '\t',quote = F,row.names = F)

# Create a data frame with all 3 possible comparisons between types of subunits for SEPV gnomAD v3
results=data.frame()
results=rbind(results,wilcox_test_func(SEPV.gnomAD.v3.Autosomal,SEPV.gnomAD.v3.Mitochondrial))
results=rbind(results,wilcox_test_func(SEPV.gnomAD.v3.Autosomal,SEPV.gnomAD.v3.X))
results=rbind(results,wilcox_test_func(SEPV.gnomAD.v3.Mitochondrial,SEPV.gnomAD.v3.X))

# Write tsv table with the results
write.table(results,'./Output/Whole_RC_SEPV_analysis/SEPV.gnomAD.v3.between.subunit.type.tsv',sep = '\t',quote = F,row.names = F)


###########################################################################
#### Compare SEPV (both meanings)  between types of subunits within CI ####
###########################################################################

# SEPV (both SEPV.gnomAD.v3 & SEPV.msa.interspecies) in autosomal in CI
SEPV.gnomAD.v3.CI.Autosomal=c(df$SEPV.gnomAD.v3[which(df$Chromosome=='Autosomal' & df$RC=='CI')],df$SEPV.gnomAD.v3[df$Subunit=='ENSG00000004779'])
SEPV.msa.interspecies.CI.Autosomal=c(df$SEPV.msa.interspecies[which(df$Chromosome=='Autosomal' & df$RC=='CI')],df$SEPV.msa.interspecies[df$Subunit=='ENSG00000004779'])

# SEPV (both SEPV.gnomAD.v3 & SEPV.msa.interspecies) in Mitochondrial in CI
SEPV.gnomAD.v3.CI.Mitochondrial=df$SEPV.gnomAD.v3[which(df$Chromosome=='Mitochondrial' & df$RC=='CI')]
SEPV.msa.interspecies.CI.Mitochondrial=df$SEPV.msa.interspecies[which(df$Chromosome=='Mitochondrial' & df$RC=='CI')]

# SEPV (both SEPV.gnomAD.v3 & SEPV.msa.interspecies) in X-linked in CI
SEPV.gnomAD.v3.CI.X=df$SEPV.gnomAD.v3[which(df$Chromosome=='X' & df$RC=='CI')]
SEPV.msa.interspecies.CI.X=df$SEPV.msa.interspecies[which(df$Chromosome=='X' & df$RC=='CI')]

# Create a data frame with all 3 possible comparisons between types of subunits for SEPV msa interspecies in CI
results=data.frame()
results=rbind(results,wilcox_test_func(SEPV.msa.interspecies.CI.Autosomal,SEPV.msa.interspecies.CI.Mitochondrial))
results=rbind(results,wilcox_test_func(SEPV.msa.interspecies.CI.Autosomal,SEPV.msa.interspecies.CI.X))
results=rbind(results,wilcox_test_func(SEPV.msa.interspecies.CI.Mitochondrial,SEPV.msa.interspecies.CI.X))

# Write tsv table with the results
write.table(results,'./Output/Whole_RC_SEPV_analysis/SEPV.msa.interspecies.between.subunit.type.in.CI.tsv',sep = '\t',quote = F,row.names = F)

# Create a data frame with all 3 possible comparisons between types of subunits for SEPV gnomAD v3 in CI
results=data.frame()
results=rbind(results,wilcox_test_func(SEPV.gnomAD.v3.CI.Autosomal,SEPV.gnomAD.v3.CI.Mitochondrial))
results=rbind(results,wilcox_test_func(SEPV.gnomAD.v3.CI.Autosomal,SEPV.gnomAD.v3.CI.X))
results=rbind(results,wilcox_test_func(SEPV.gnomAD.v3.CI.Mitochondrial,SEPV.gnomAD.v3.CI.X))

# Write tsv table with the results
write.table(results,'./Output/Whole_RC_SEPV_analysis/SEPV.gnomAD.v3.between.subunit.type.in.CI.tsv',sep = '\t',quote = F,row.names = F)


###################################################################
#### Compare SEPV (both meanings) by types of subunits in CIII ####
###################################################################

# SEPV (both SEPV.gnomAD.v3 & SEPV.msa.interspecies) in autosomal in CIII
SEPV.gnomAD.v3.CIII.Autosomal=c(df$SEPV.gnomAD.v3[which(df$Chromosome=='Autosomal' & df$RC=='CIII')],df$SEPV.gnomAD.v3[which(df$Chromosome=='Autosomal' & df$RC=='CIII')])
SEPV.msa.interspecies.CIII.Autosomal=c(df$SEPV.msa.interspecies[which(df$Chromosome=='Autosomal' & df$RC=='CIII')],df$SEPV.msa.interspecies[which(df$Chromosome=='Autosomal' & df$RC=='CIII')])

# SEPV (both SEPV.gnomAD.v3 & SEPV.msa.interspecies) in Mitochondrial in CIII
SEPV.gnomAD.v3.CIII.Mitochondrial=c(df$SEPV.gnomAD.v3[which(df$Chromosome=='Mitochondrial' & df$RC=='CIII')],df$SEPV.gnomAD.v3[which(df$Chromosome=='Mitochondrial' & df$RC=='CIII')])
SEPV.msa.interspecies.CIII.Mitochondrial=c(df$SEPV.msa.interspecies[which(df$Chromosome=='Mitochondrial' & df$RC=='CIII')],df$SEPV.gnomAD.v3[which(df$Chromosome=='Mitochondrial' & df$RC=='CIII')])

# IMPORTANT: In CIII there are not X-linked subunits

# Create a data frame with all 3 possible comparisons between types of subunits for SEPV msa interspecies
results=data.frame()
results=rbind(results,wilcox_test_func(SEPV.msa.interspecies.CIII.Autosomal,SEPV.msa.interspecies.CIII.Mitochondrial))

# Write tsv table with the results
write.table(results,'./Output/Whole_RC_SEPV_analysis/SEPV.gnomAD.v3.between.subunit.type.in.CIII.tsv',sep = '\t',quote = F,row.names = F)

# Create a data frame with all 3 possible comparisons between types of subunits for SEPV gnomAD v3 in CIII
results=data.frame()
results=rbind(results,wilcox_test_func(SEPV.gnomAD.v3.CIII.Autosomal,SEPV.gnomAD.v3.CIII.Mitochondrial))

# Write tsv table with the results
write.table(results,'./Output/Whole_RC_SEPV_analysis/SEPV.msa.interspecies.between.subunit.type.in.CIII.tsv',sep = '\t',quote = F,row.names = F)


##################################################################
#### Compare SEPV (both meanings) by types of subunits in CIV ####
##################################################################

# SEPV (both SEPV.gnomAD.v3 & SEPV.msa.interspecies) in autosomal in CIV
SEPV.gnomAD.v3.CIV.Autosomal=df$SEPV.gnomAD.v3[which(df$Chromosome=='Autosomal' & df$RC=='CIV')]
SEPV.msa.interspecies.CIV.Autosomal=df$SEPV.msa.interspecies[which(df$Chromosome=='Autosomal' & df$RC=='CIV')]

# SEPV (both SEPV.gnomAD.v3 & SEPV.msa.interspecies) in Mitochondrial in CIV
SEPV.gnomAD.v3.CIV.Mitochondrial=df$SEPV.gnomAD.v3[which(df$Chromosome=='Mitochondrial' & df$RC=='CIV')]
SEPV.msa.interspecies.CIV.Mitochondrial=df$SEPV.msa.interspecies[which(df$Chromosome=='Mitochondrial' & df$RC=='CIV')]

# SEPV (both SEPV.gnomAD.v3 & SEPV.msa.interspecies) in X-linked in CIV
SEPV.gnomAD.v3.CIV.X=df$SEPV.gnomAD.v3[which(df$Chromosome=='X' & df$RC=='CIV')]
SEPV.msa.interspecies.CIV.X=df$SEPV.msa.interspecies[which(df$Chromosome=='X' & df$RC=='CIV')]

# Create a data frame with all 3 possible comparisons between types of subunits for SEPV msa interspecies
results=data.frame()
results=rbind(results,wilcox_test_func(SEPV.msa.interspecies.CIV.Autosomal,SEPV.msa.interspecies.CIV.Mitochondrial))
results=rbind(results,wilcox_test_func(SEPV.msa.interspecies.CIV.Autosomal,SEPV.msa.interspecies.CIV.X))
results=rbind(results,wilcox_test_func(SEPV.msa.interspecies.CIV.Mitochondrial,SEPV.msa.interspecies.CIV.X))

# Write tsv table with the results
write.table(results,'./Output/Whole_RC_SEPV_analysis/SEPV.msa.interspecies.between.subunit.type.in.CIV.tsv',sep = '\t',quote = F,row.names = F)

# Create a data frame with all 3 possible comparisons between types of subunits for SEPV gnomAD v3 in CIV
results=data.frame()
results=rbind(results,wilcox_test_func(SEPV.gnomAD.v3.CIV.Autosomal,SEPV.gnomAD.v3.CIV.Mitochondrial))
results=rbind(results,wilcox_test_func(SEPV.gnomAD.v3.CIV.Autosomal,SEPV.gnomAD.v3.CIV.X))
results=rbind(results,wilcox_test_func(SEPV.gnomAD.v3.CIV.Mitochondrial,SEPV.gnomAD.v3.CIV.X))

# Write tsv table with the results
write.table(results,'./Output/Whole_RC_SEPV_analysis/SEPV.gnomAD.v3.between.subunit.type.in.CIV.tsv',sep = '\t',quote = F,row.names = F)


#################################################################
#### Compare SEPV (both meanings) in types of subunits in CV ####
#################################################################

# SEPV (both SEPV.gnomAD.v3 & SEPV.msa.interspecies) in autosomal in CV
SEPV.gnomAD.v3.CV.Autosomal=c(df$SEPV.gnomAD.v3[which(df$Chromosome=='Autosomal' & df$RC=='CV')],rep(df$SEPV.gnomAD.v3[df$Subunit=='ENSG00000152234'],2), 
                              rep(df$SEPV.gnomAD.v3[df$Subunit=='ENSG00000110955'],2), rep(df$SEPV.gnomAD.v3[df$Subunit=='ENSG00000159199'],7),
                              rep(df$SEPV.gnomAD.v3[df$Subunit=='ENSG00000135390'],7),rep(df$SEPV.gnomAD.v3[df$Subunit=='ENSG00000135390'],7))
SEPV.msa.interspecies.CV.Autosomal=c(df$SEPV.msa.interspecies[which(df$Chromosome=='Autosomal' & df$RC=='CV')],rep(df$SEPV.msa.interspecies[df$Subunit=='ENSG00000152234'],2), 
                                     rep(df$SEPV.msa.interspecies[df$Subunit=='ENSG00000110955'],2), rep(df$SEPV.msa.interspecies[df$Subunit=='ENSG00000159199'],7),
                                     rep(df$SEPV.msa.interspecies[df$Subunit=='ENSG00000135390'],7),rep(df$SEPV.msa.interspecies[df$Subunit=='ENSG00000135390'],7))

# SEPV (both SEPV.gnomAD.v3 & SEPV.msa.interspecies) in Mitochondrial in CV
SEPV.gnomAD.v3.CV.Mitochondrial=df$SEPV.gnomAD.v3[which(df$Chromosome=='Mitochondrial' & df$RC=='CV')]
SEPV.msa.interspecies.CV.Mitochondrial=df$SEPV.msa.interspecies[which(df$Chromosome=='Mitochondrial' & df$RC=='CV')]


# IMPORTANT: In CV there are no X-linked subunits

# Create a data frame with the only possible comparison between types of subunits for SEPV msa interspecies
results=data.frame()
results=rbind(results,wilcox_test_func(SEPV.msa.interspecies.CV.Autosomal,SEPV.msa.interspecies.CV.Mitochondrial))

# Write tsv table with the results
write.table(results,'./Output/Whole_RC_SEPV_analysis/SEPV.msa.interspecies.between.subunit.type.in.CV.tsv',sep = '\t',quote = F,row.names = F)

# Create a data frame with the only possible comparison between types of subunits for SEPV gnomAD v3 in CV
results=data.frame()
results=rbind(results,wilcox_test_func(SEPV.gnomAD.v3.CV.Autosomal,SEPV.gnomAD.v3.CV.Mitochondrial))

# Write tsv table with the results
write.table(results,'./Output/Whole_RC_SEPV_analysis/SEPV.gnomAD.v3.between.subunit.type.in.CV.tsv',sep = '\t',quote = F,row.names = F)


########################################################################
#### Compare SEPV (both meanings) between RCs in autosomal subunits ####
########################################################################

# SEPV (both SEPV.gnomAD.v3 & SEPV.msa.interspecies) CI in autosomal subunits
CI.SEPV.gnomAD.v3.in.autosomal=c(df$SEPV.gnomAD.v3[df$RC=='CI' & df$Chromosome=='Autosomal'],df$SEPV.gnomAD.v3[df$Subunit=='ENSG00000004779'])
CI.SEPV.msa.interspecies.in.autosomal=c(df$SEPV.msa.interspecies[df$RC=='CI' & df$Chromosome=='Autosomal'],df$SEPV.msa.interspecies[df$Subunit=='ENSG00000004779'])

# SEPV (both SEPV.gnomAD.v3 & SEPV.msa.interspecies) CIII in autosomal subunits
CII.SEPV.gnomAD.v3.in.autosomal=df$SEPV.gnomAD.v3[df$RC=='CII' & df$Chromosome=='Autosomal']
CII.SEPV.msa.interspecies.in.autosomal=df$SEPV.msa.interspecies[df$RC=='CII'& df$Chromosome=='Autosomal']

# SEPV (both SEPV.gnomAD.v3 & SEPV.msa.interspecies) CIII2I in autosomal subunits
CIII.SEPV.gnomAD.v3.in.autosomal=c(df$SEPV.gnomAD.v3[df$RC=='CIII' & df$Chromosome=='Autosomal'],df$SEPV.gnomAD.v3[df$RC=='CIII' & df$Chromosome=='Autosomal'])
CIII.SEPV.msa.interspecies.in.autosomal=c(df$SEPV.msa.interspecies[df$RC=='CIII'& df$Chromosome=='Autosomal'],df$SEPV.msa.interspecies[df$RC=='CIII'& df$Chromosome=='Autosomal'])

# SEPV (both SEPV.gnomAD.v3 & SEPV.msa.interspecies) CIVI in autosomal subunits
CIV.SEPV.gnomAD.v3.in.autosomal=df$SEPV.gnomAD.v3[df$RC=='CIV' & df$Chromosome=='Autosomal']
CIV.SEPV.msa.interspecies.in.autosomal=df$SEPV.msa.interspecies[df$RC=='CIV' & df$Chromosome=='Autosomal']


# SEPV (both SEPV.gnomAD.v3 & SEPV.msa.interspecies) CVI in autosomal subunits
CV.SEPV.gnomAD.v3.in.autosomal=c(df$SEPV.gnomAD.v3[df$RC=='CV' & df$Chromosome=='Autosomal'],rep(df$SEPV.gnomAD.v3[df$Subunit=='ENSG00000152234'],2),
                                 rep(df$SEPV.gnomAD.v3[df$Subunit=='ENSG00000110955'],2), rep(df$SEPV.gnomAD.v3[df$Subunit=='ENSG00000159199'],7), 
                                 rep(df$SEPV.gnomAD.v3[df$Subunit=='ENSG00000135390'],7), rep(df$SEPV.gnomAD.v3[df$Subunit=='ENSG00000154518'],7))
CV.SEPV.msa.interspecies.in.autosomal=c(df$SEPV.msa.interspecies[df$RC=='CV' & df$Chromosome=='Autosomal'],rep(df$SEPV.msa.interspecies[df$Subunit=='ENSG00000152234'],2),
                                        rep(df$SEPV.msa.interspecies[df$Subunit=='ENSG00000110955'],2), rep(df$SEPV.msa.interspecies[df$Subunit=='ENSG00000159199'],7), 
                                        rep(df$SEPV.msa.interspecies[df$Subunit=='ENSG00000135390'],7), rep(df$SEPV.msa.interspecies[df$Subunit=='ENSG00000154518'],7))


# Create a data frame with all 10 possible comparisons between RCs for SEPV msa interspecies for autosomal subunits
results=data.frame()
results=rbind(results,wilcox_test_func(CI.SEPV.msa.interspecies.in.autosomal,CII.SEPV.msa.interspecies.in.autosomal))
results=rbind(results,wilcox_test_func(CI.SEPV.msa.interspecies.in.autosomal,CIII.SEPV.msa.interspecies.in.autosomal))
results=rbind(results,wilcox_test_func(CI.SEPV.msa.interspecies.in.autosomal,CIV.SEPV.msa.interspecies.in.autosomal))
results=rbind(results,wilcox_test_func(CI.SEPV.msa.interspecies.in.autosomal,CV.SEPV.msa.interspecies.in.autosomal))
results=rbind(results,wilcox_test_func(CII.SEPV.msa.interspecies.in.autosomal,CIII.SEPV.msa.interspecies.in.autosomal))
results=rbind(results,wilcox_test_func(CII.SEPV.msa.interspecies.in.autosomal,CIV.SEPV.msa.interspecies.in.autosomal))
results=rbind(results,wilcox_test_func(CII.SEPV.msa.interspecies.in.autosomal,CV.SEPV.msa.interspecies.in.autosomal))
results=rbind(results,wilcox_test_func(CIII.SEPV.msa.interspecies.in.autosomal,CIV.SEPV.msa.interspecies.in.autosomal))
results=rbind(results,wilcox_test_func(CIII.SEPV.msa.interspecies.in.autosomal,CV.SEPV.msa.interspecies.in.autosomal))
results=rbind(results,wilcox_test_func(CIV.SEPV.msa.interspecies.in.autosomal,CV.SEPV.msa.interspecies.in.autosomal))

# Write tsv table with the results
write.table(results,'./Output/Whole_RC_SEPV_analysis/SEPV.msa.interspecies.in.autosomal.between.RCs.tsv',sep = '\t',quote = F,row.names = F)

# Create a data frame with all 10 possible comparisons between RCs for SEPV msa interspecies for autosomal subunits
results=data.frame()
results=rbind(results,wilcox_test_func(CI.SEPV.gnomAD.v3.in.autosomal,CII.SEPV.gnomAD.v3.in.autosomal))
results=rbind(results,wilcox_test_func(CI.SEPV.gnomAD.v3.in.autosomal,CIII.SEPV.gnomAD.v3.in.autosomal))
results=rbind(results,wilcox_test_func(CI.SEPV.gnomAD.v3.in.autosomal,CIV.SEPV.gnomAD.v3.in.autosomal))
results=rbind(results,wilcox_test_func(CI.SEPV.gnomAD.v3.in.autosomal,CV.SEPV.gnomAD.v3.in.autosomal))
results=rbind(results,wilcox_test_func(CII.SEPV.gnomAD.v3.in.autosomal,CIII.SEPV.gnomAD.v3.in.autosomal))
results=rbind(results,wilcox_test_func(CII.SEPV.gnomAD.v3.in.autosomal,CIV.SEPV.gnomAD.v3.in.autosomal))
results=rbind(results,wilcox_test_func(CII.SEPV.gnomAD.v3.in.autosomal,CV.SEPV.gnomAD.v3.in.autosomal))
results=rbind(results,wilcox_test_func(CIII.SEPV.gnomAD.v3.in.autosomal,CIV.SEPV.gnomAD.v3.in.autosomal))
results=rbind(results,wilcox_test_func(CIII.SEPV.gnomAD.v3.in.autosomal,CV.SEPV.gnomAD.v3.in.autosomal))
results=rbind(results,wilcox_test_func(CIV.SEPV.gnomAD.v3.in.autosomal,CV.SEPV.gnomAD.v3.in.autosomal))

# Write tsv table with the results
write.table(results,'./Output/Whole_RC_SEPV_analysis/SEPV.gnomAD.v3.in.autosomal.between.RCs.tsv',sep = '\t',quote = F,row.names = F)


###############################################
#### SEPV by RCs in mitochondrial subunits ####
###############################################

# SEPV (both SEPV.gnomAD.v3 & SEPV.msa.interspecies) CI I in mitichondrial subunits
CI.SEPV.gnomAD.v3.in.Mitochondrial=df$SEPV.gnomAD.v3[df$RC=='CI' & df$Chromosome=='Mitochondrial']
CI.SEPV.msa.interspecies.in.Mitochondrial=df$SEPV.msa.interspecies[df$RC=='CI' & df$Chromosome=='Mitochondrial']


# SEPV (both SEPV.gnomAD.v3 & SEPV.msa.interspecies) CIII2 in mitichondrial subunits
CIII.SEPV.gnomAD.v3.in.Mitochondrial=c(df$SEPV.gnomAD.v3[df$RC=='CIII' & df$Chromosome=='Mitochondrial'],df$SEPV.gnomAD.v3[df$RC=='CIII' & df$Chromosome=='Mitochondrial'])
CIII.SEPV.msa.interspecies.in.Mitochondrial=c(df$SEPV.msa.interspecies[df$RC=='CIII'& df$Chromosome=='Mitochondrial'],df$SEPV.msa.interspecies[df$RC=='CIII'& df$Chromosome=='Mitochondrial'])

# SEPV (both SEPV.gnomAD.v3 & SEPV.msa.interspecies) CIV in mitichondrial subunits
CIV.SEPV.gnomAD.v3.in.Mitochondrial=df$SEPV.gnomAD.v3[df$RC=='CIV' & df$Chromosome=='Mitochondrial']
CIV.SEPV.msa.interspecies.in.Mitochondrial=df$SEPV.msa.interspecies[df$RC=='CIV' & df$Chromosome=='Mitochondrial']


# SEPV (both SEPV.gnomAD.v3 & SEPV.msa.interspecies) CV in mitichondrial subunits
CV.SEPV.gnomAD.v3.in.Mitochondrial=df$SEPV.gnomAD.v3[df$RC=='CV' & df$Chromosome=='Mitochondrial']
CV.SEPV.msa.interspecies.in.Mitochondrial=df$SEPV.msa.interspecies[df$RC=='CV' & df$Chromosome=='Mitochondrial']



# Create a data frame with all 6 possible comparisons between RCs for SEPV msa interspecies for mitochondrial subunits
results=data.frame()
results=rbind(results,wilcox_test_func(CI.SEPV.msa.interspecies.in.Mitochondrial,CIII.SEPV.msa.interspecies.in.Mitochondrial))
results=rbind(results,wilcox_test_func(CI.SEPV.msa.interspecies.in.Mitochondrial,CIV.SEPV.msa.interspecies.in.Mitochondrial))
results=rbind(results,wilcox_test_func(CI.SEPV.msa.interspecies.in.Mitochondrial,CV.SEPV.msa.interspecies.in.Mitochondrial))
results=rbind(results,wilcox_test_func(CIII.SEPV.msa.interspecies.in.Mitochondrial,CIV.SEPV.msa.interspecies.in.Mitochondrial))
results=rbind(results,wilcox_test_func(CIII.SEPV.msa.interspecies.in.Mitochondrial,CV.SEPV.msa.interspecies.in.Mitochondrial))
results=rbind(results,wilcox_test_func(CIV.SEPV.msa.interspecies.in.Mitochondrial,CV.SEPV.msa.interspecies.in.Mitochondrial))

# Write tsv table with the results
write.table(results,'./Output/Whole_RC_SEPV_analysis/SEPV.msa.interspecies.in.Mitochondrial.between.RCs.tsv',sep = '\t',quote = F,row.names = F)

# Create a data frame with all 6 possible comparisons between RCs for SEPV msa interspecies for Mitochondrial subunits
results=data.frame()
results=rbind(results,wilcox_test_func(CI.SEPV.gnomAD.v3.in.Mitochondrial,CIII.SEPV.gnomAD.v3.in.Mitochondrial))
results=rbind(results,wilcox_test_func(CI.SEPV.gnomAD.v3.in.Mitochondrial,CIV.SEPV.gnomAD.v3.in.Mitochondrial))
results=rbind(results,wilcox_test_func(CI.SEPV.gnomAD.v3.in.Mitochondrial,CV.SEPV.gnomAD.v3.in.Mitochondrial))
results=rbind(results,wilcox_test_func(CIII.SEPV.gnomAD.v3.in.Mitochondrial,CIV.SEPV.gnomAD.v3.in.Mitochondrial))
results=rbind(results,wilcox_test_func(CIII.SEPV.gnomAD.v3.in.Mitochondrial,CV.SEPV.gnomAD.v3.in.Mitochondrial))
results=rbind(results,wilcox_test_func(CIV.SEPV.gnomAD.v3.in.Mitochondrial,CV.SEPV.gnomAD.v3.in.Mitochondrial))

# Write tsv table with the results
write.table(results,'./Output/Whole_RC_SEPV_analysis/SEPV.gnomAD.v3.in.Mitochondrial.between.RCs.tsv',sep = '\t',quote = F,row.names = F)


###################################
#### SEPV by RCs in X subunits ####
###################################

# SEPV (both SEPV.gnomAD.v3 & SEPV.msa.interspecies) CI 
CI.SEPV.gnomAD.v3.in.X=df$SEPV.gnomAD.v3[df$RC=='CI' & df$Chromosome=='X']
CI.SEPV.msa.interspecies.in.X=df$SEPV.msa.interspecies[df$RC=='CI' & df$Chromosome=='X']

# SEPV (both SEPV.gnomAD.v3 & SEPV.msa.interspecies) CIV
CIV.SEPV.gnomAD.v3.in.X=df$SEPV.gnomAD.v3[df$RC=='CIV' & df$Chromosome=='X']
CIV.SEPV.msa.interspecies.in.X=df$SEPV.msa.interspecies[df$RC=='CIV' & df$Chromosome=='X']

# Create a data frame with unique possible comparison between RCs for SEPV msa interspecies for X subunits
results=data.frame()
results=rbind(results,wilcox_test_func(CI.SEPV.msa.interspecies.in.X,CIV.SEPV.msa.interspecies.in.X))

# Write tsv table with the results
write.table(results,'./Output/Whole_RC_SEPV_analysis/SEPV.msa.interspecies.in.X.between.RCs.tsv',sep = '\t',quote = F,row.names = F)


# Create a data frame with unique possible comparisons between RCs for SEPV msa interspecies for X subunits
results=data.frame()
results=rbind(results,wilcox_test_func(CI.SEPV.gnomAD.v3.in.X,CIV.SEPV.gnomAD.v3.in.X))

# Write tsv table with the results
write.table(results,'./Output/Whole_RC_SEPV_analysis/SEPV.gnomAD.v3.in.X.between.RCs.tsv',sep = '\t',quote = F,row.names = F)


##################################################
### GETS IFR FROM IN SILICO MSTRUCTURES OF RCs ###
##################################################

# Define class for object S4 where all information will be stored
setClass("Interfaces", representation(name = "vector", Ens.subunit.query = "vector",Ens.subunit.target='vector',
                                      subunit.query = "vector",subunit.target='vector',
                                      Interface.residues.query='list',genome.query="vector",
                                      genome.target="vector",SEPV.gnomAD.v3='list',SEPV.msa.interspecies='list'))

# Create object S4 of class interfaces
p1=new("Interfaces",name=character(),Ens.subunit.query=character(),
       Ens.subunit.target=character(),subunit.query = character(),subunit.target=character(),
       Interface.residues.query=list(),SEPV.msa.interspecies=list(),SEPV.gnomAD.v3=list(),
       genome.query=character(),genome.target=character())



# List of structure files and their respective tsv files that allow mapping each chain with its respective gene
path='./input'
l1=list.files(path,pattern = '.Chain.gene.key.tsv')
l2=list.files(path,pattern = '.pdb')

# Iterate files to be analyzed 
for (j in 1:length(l1)){
  
  # Load pdb structure and its respective tsv file to map chain Vs gene
  key=read.delim(paste0(path,'/',l1[j]),sep = '\t',stringsAsFactors = F)
  pdb <- read.pdb(paste0(path,'/',l2[j]))
  
  # Get name of studied respiratory complex (RC)
  name=gsub('.pdb','',l2[j])
  
  # Define all possible combinations of pairs chains (query-target) in the structure to be evaluated in search for IFRs
  grid=expand.grid(unique(pdb$atom$chain),unique(pdb$atom$chain))
  grid=grid[which(grid$Var1!=grid$Var2),]
  
  # Annotate S4 object with RC information
  p1@name=c(p1@name,name)
  
  # Annotate the S4 object with the Ensembl ID and the symbol of the gene representing the evaluated chains
  p1@Ens.subunit.query=c(p1@Ens.subunit.query,sapply(grid[,1],function(x)return(key[key[["Chain"]]==x,"ENSEMBL_ID"])))
  p1@Ens.subunit.target=c(p1@Ens.subunit.target,sapply(grid[,2],function(x)return(key[key[["Chain"]]==x,"ENSEMBL_ID"])))
  p1@subunit.query=c(p1@subunit.query,sapply(grid[,1],function(x)return(key[key[["Chain"]]==x,"Gene"])))
  p1@subunit.target=c(p1@subunit.target,sapply(grid[,2],function(x)return(key[key[["Chain"]]==x,"Gene"])))
  
  # Annotate the IFRs between chains with their coresponding values for SEPV.gnomAD.v3 and SEPV.msa.interspecies
  p1@Interface.residues.query=c(p1@Interface.residues.query,mclapply(1:nrow(grid),function(i){binding.site(pdb, a.inds=atom.select(pdb, chain=as.character(grid[i,1])), b.inds=atom.select(pdb, chain=as.character(grid[i,2])), verbose = F)$resno},mc.cores = 11))
  
  p1@SEPV.gnomAD.v3=c(p1@SEPV.gnomAD.v3,mclapply(1:nrow(grid),function(i){
    res=binding.site(pdb, a.inds=atom.select(pdb, chain=as.character(grid[i,1])), b.inds=atom.select(pdb, chain=as.character(grid[i,2])), verbose = F)$resno
    gene=key$ENSEMBL_ID[which(key[["Chain"]]==grid[i,1])]
    if (length(res)>0){
      df.IFR=df %>% filter(Subunit==gene & Prot.pos %in% res) %>% arrange(match(Prot.pos,res))
      v2=df.IFR[["SEPV.gnomAD.v3"]]
    }
    else{v2=NA}
    return(v2)
  },mc.cores = 11))
  
  p1@SEPV.msa.interspecies=c(p1@SEPV.msa.interspecies,mclapply(1:nrow(grid),function(i){
    res=binding.site(pdb, a.inds=atom.select(pdb, chain=as.character(grid[i,1])), b.inds=atom.select(pdb, chain=as.character(grid[i,2])), verbose = F)$resno
    gene=key$ENSEMBL_ID[which(key[["Chain"]]==grid[i,1])]
    if (length(res)>0){
      df.IFR=df %>% filter(Subunit==gene & Prot.pos %in% res) %>% arrange(match(Prot.pos,res))
      v2=df.IFR[["SEPV.msa.interspecies"]]
    }
    else{v2=NA}
    return(v2)
  },mc.cores = 11))
  
}

# Annotate IFRs between chains with their values of encoding chromosomes for involved chains at the binding site
p1@genome.query=p1@Ens.subunit.query
p1@genome.query[p1@genome.query %in% X.genes]='X'
p1@genome.query[p1@genome.query %in% Mt.genes]='Mitochondrial'
p1@genome.query[p1@genome.query %in% c("X","Mitochondrial") == F]='Autosomal'

p1@genome.target=p1@Ens.subunit.target
p1@genome.target[p1@genome.target %in% X.genes]='X'
p1@genome.target[p1@genome.target %in% Mt.genes]='Mitochondrial'
p1@genome.target[p1@genome.target %in% c("X","Mitochondrial") == F]='Autosomal'

# Eliminates all non-neighboring combinations between chains 
cond=unlist(lapply(p1@Interface.residues.query,function(x){return(length(x))}))>0
p1@Ens.subunit.query=p1@Ens.subunit.query[cond]
p1@Ens.subunit.target=p1@Ens.subunit.target[cond]
p1@subunit.query=p1@subunit.query[cond]
p1@subunit.target=p1@subunit.target[cond]
p1@genome.query=p1@genome.query[cond]
p1@genome.target=p1@genome.target[cond]
p1@Interface.residues.query=p1@Interface.residues.query[cond]
p1@SEPV.gnomAD.v3=p1@SEPV.gnomAD.v3[cond]
p1@SEPV.msa.interspecies=p1@SEPV.msa.interspecies[cond]

# Save S4 object containing information for IFRs
saveRDS(p1,"./Output/IFRs_SEPV_analysis/IFRs.rds")

##################################################################################################
# Comparison of conservation and population variability meanings of SEPV between IFRs Vs NO.IFRs #
##################################################################################################

# Get SEPV.msa.interspecies values for IFRs from in silico models
IFRs.SEPV.msa.interspecies=unlist(p1@SEPV.msa.interspecies)
names(IFRs.SEPV.msa.interspecies)=unlist(sapply(1:length(p1@SEPV.gnomAD.v3),function(x){paste0(p1@subunit.query[x],':',unlist(p1@Interface.residues.query[x]))}))
IFRs.SEPV.msa.interspecies

# Get SEPV.gnomAD.v3 values for IFRs from in silico models
IFRs.SEPV.gnomAD.v3=unlist(p1@SEPV.gnomAD.v3)
names(IFRs.SEPV.gnomAD.v3)=unlist(sapply(1:length(p1@SEPV.gnomAD.v3),function(x){paste0(p1@subunit.query[x],':',unlist(p1@Interface.residues.query[x]))}))
IFRs.SEPV.gnomAD.v3

# Get SEPV.msa.interspecies values for IFRs from in silico models for autosomal encoded genes
IFRs.SEPV.msa.interspecies.autosomal=unlist(p1@SEPV.msa.interspecies[which(p1@genome.query=='Autosomal')])
names(IFRs.SEPV.msa.interspecies.autosomal)=unlist(sapply(which(p1@genome.query=='Autosomal'),function(x){paste0(p1@subunit.query[x],':',unlist(p1@Interface.residues.query[x]))}))
IFRs.SEPV.msa.interspecies.autosomal

# Get SEPV.gnomAD.v3 values for IFRs from in silico models for autosomal encoded genes
IFRs.SEPV.gnomAD.v3.autosomal=unlist(p1@SEPV.gnomAD.v3[which(p1@genome.query=='Autosomal')])
names(IFRs.SEPV.gnomAD.v3.autosomal)=unlist(sapply(which(p1@genome.query=='Autosomal'),function(x){paste0(p1@subunit.query[x],':',unlist(p1@Interface.residues.query[x]))}))
IFRs.SEPV.gnomAD.v3.autosomal

# Get SEPV.msa.interspecies values for IFRs from in silico models for Mitochondrial encoded genes
IFRs.SEPV.msa.interspecies.Mitochondrial=unlist(p1@SEPV.msa.interspecies[which(p1@genome.query=='Mitochondrial')])
names(IFRs.SEPV.msa.interspecies.Mitochondrial)=unlist(sapply(which(p1@genome.query=='Mitochondrial'),function(x){paste0(p1@subunit.query[x],':',unlist(p1@Interface.residues.query[x]))}))
IFRs.SEPV.msa.interspecies.Mitochondrial

# Get SEPV.gnomAD.v3 values for IFRs from in silico models for Mitochondrial encoded genes
IFRs.SEPV.gnomAD.v3.Mitochondrial=unlist(p1@SEPV.gnomAD.v3[which(p1@genome.query=='Mitochondrial')])
names(IFRs.SEPV.gnomAD.v3.Mitochondrial)=unlist(sapply(which(p1@genome.query=='Mitochondrial'),function(x){paste0(p1@subunit.query[x],':',unlist(p1@Interface.residues.query[x]))}))
IFRs.SEPV.gnomAD.v3.Mitochondrial

# Get SEPV.msa.interspecies values for IFRs from in silico models for X-linked genes
IFRs.SEPV.msa.interspecies.X=unlist(p1@SEPV.msa.interspecies[which(p1@genome.query=='X')])
names(IFRs.SEPV.msa.interspecies.X)=unlist(sapply(which(p1@genome.query=='X'),function(x){paste0(p1@subunit.query[x],':',unlist(p1@Interface.residues.query[x]))}))
IFRs.SEPV.msa.interspecies.X

# Get SEPV.gnomAD.v3 values for IFRs from in silico models for X-linked genes
IFRs.SEPV.gnomAD.v3.X=unlist(p1@SEPV.gnomAD.v3[which(p1@genome.query=='X')])
names(IFRs.SEPV.gnomAD.v3.X)=unlist(sapply(which(p1@genome.query=='X'),function(x){paste0(p1@subunit.query[x],':',unlist(p1@Interface.residues.query[x]))}))
IFRs.SEPV.gnomAD.v3.X


# Get Non.IFRs considering their stoichiometry in our in silico models
dfni=df[df$Subunit %in% unique(p1@Ens.subunit.query),]
dfni=rbind(dfni,dfni[dfni$Subunit=='ENSG00000004779',], dfni[rep(rownames(dfni[dfni$Subunit=='ENSG00000152234',]),2),], dfni[dfni$RC=='CIII',],
         dfni[rep(rownames(dfni[dfni$Subunit=='ENSG00000110955',]),2),], dfni[rep(rownames(dfni[dfni$Subunit=='ENSG00000159199',]),7),])

# Get SEPV.gnomAD.v3 values for IFRs from in silico models
NO.IFRs.SEPV.msa.interspecies=dfni$SEPV.msa.interspecies[paste0(dfni$Subunit,':',dfni$Prot.pos) %in% names(IFRs.SEPV.msa.interspecies)==F]

# Get SEPV.msa.interspecies values for IFRs from in silico models
NO.IFRs.SEPV.gnomAD.v3=dfni$SEPV.gnomAD.v3[paste0(dfni$Subunit,':',dfni$Prot.pos) %in% names(IFRs.SEPV.msa.interspecies)==F]

# Get SEPV.gnomAD.v3 values for IFRs from in silico models for autosomal genes
NO.IFRs.SEPV.gnomAD.v3.autosomal=dfni$SEPV.gnomAD.v3[paste0(dfni$Subunit,':',dfni$Prot.pos) %in% names(IFRs.SEPV.gnomAD.v3.autosomal)==F & dfni$Subunit %in% p1@Ens.subunit.query[p1@genome.query=='Autosomal']]

# Get SEPV.msa.interspecies values for IFRs from in silico models for autosomal genes
NO.IFRs.SEPV.msa.interspecies.autosomal=dfni$SEPV.msa.interspecies[paste0(dfni$Subunit,':',dfni$Prot.pos) %in% names(IFRs.SEPV.msa.interspecies.autosomal)==F & dfni$Subunit %in% p1@Ens.subunit.query[p1@genome.query=='Autosomal']]

# Get SEPV.gnomAD.v3 values for IFRs from in silico models for mitochondrial genes
NO.IFRs.SEPV.gnomAD.v3.Mitochondrial=dfni$SEPV.gnomAD.v3[paste0(dfni$Subunit,':',dfni$Prot.pos) %in% names(IFRs.SEPV.gnomAD.v3.Mitochondrial)==F & dfni$Subunit %in% p1@Ens.subunit.query[p1@genome.query=='Mitochondrial']]

# Get SEPV.msa.interspecies values for IFRs from in silico models for mitochondrial genes
NO.IFRs.SEPV.msa.interspecies.Mitochondrial=dfni$SEPV.msa.interspecies[paste0(dfni$Subunit,':',dfni$Prot.pos) %in% names(IFRs.SEPV.msa.interspecies.Mitochondrial)==F & dfni$Subunit %in% p1@Ens.subunit.query[p1@genome.query=='Mitochondrial']]

# Get SEPV.gnomAD.v3 values for IFRs from in silico models for X genes
NO.IFRs.SEPV.gnomAD.v3.X=dfni$SEPV.gnomAD.v3[paste0(dfni$Subunit,':',dfni$Prot.pos) %in% names(IFRs.SEPV.gnomAD.v3.X)==F & dfni$Subunit %in% p1@Ens.subunit.query[p1@genome.query=='X']]

# Get SEPV.msa.interspecies values for IFRs from in silico models for X genes
NO.IFRs.SEPV.msa.interspecies.X=dfni$SEPV.msa.interspecies[paste0(dfni$Subunit,':',dfni$Prot.pos) %in% names(IFRs.SEPV.msa.interspecies.X)==F & dfni$Subunit %in% p1@Ens.subunit.query[p1@genome.query=='X']]

# Create a data frame with all 10 possible comparisons IFR Vs NO-IFRs for for SEPV msa interspecies
results=data.frame()
results=rbind(results,wilcox_test_func(NO.IFRs.SEPV.msa.interspecies,IFRs.SEPV.msa.interspecies))
results=rbind(results,wilcox_test_func(NO.IFRs.SEPV.msa.interspecies.autosomal,IFRs.SEPV.msa.interspecies.autosomal))
results=rbind(results,wilcox_test_func(NO.IFRs.SEPV.msa.interspecies.Mitochondrial,IFRs.SEPV.msa.interspecies.autosomal))
results=rbind(results,wilcox_test_func(NO.IFRs.SEPV.msa.interspecies.X,IFRs.SEPV.msa.interspecies.autosomal))
results=rbind(results,wilcox_test_func(NO.IFRs.SEPV.msa.interspecies.autosomal,IFRs.SEPV.msa.interspecies.Mitochondrial))
results=rbind(results,wilcox_test_func(NO.IFRs.SEPV.msa.interspecies.Mitochondrial,IFRs.SEPV.msa.interspecies.Mitochondrial))
results=rbind(results,wilcox_test_func(NO.IFRs.SEPV.msa.interspecies.X,IFRs.SEPV.msa.interspecies.Mitochondrial))
results=rbind(results,wilcox_test_func(NO.IFRs.SEPV.msa.interspecies.autosomal,IFRs.SEPV.msa.interspecies.X))
results=rbind(results,wilcox_test_func(NO.IFRs.SEPV.msa.interspecies.Mitochondrial,IFRs.SEPV.msa.interspecies.X))
results=rbind(results,wilcox_test_func(NO.IFRs.SEPV.msa.interspecies.X,IFRs.SEPV.msa.interspecies.X))

# Write tsv table with the results
write.table(results,'./Output/IFRs_SEPV_analysis/SEPV.msa.interspecies.NO.IFRs.Vs.IFRs.by.subunit.type.tsv',sep = '\t',quote = F,row.names = F)


# Create a data frame with all 10 possible comparisons IFR Vs NO-IFRs for SEPV gnomAD v3
results=data.frame()
results=rbind(results,wilcox_test_func(NO.IFRs.SEPV.gnomAD.v3,IFRs.SEPV.gnomAD.v3))
results=rbind(results,wilcox_test_func(NO.IFRs.SEPV.gnomAD.v3.autosomal,IFRs.SEPV.gnomAD.v3.autosomal))
results=rbind(results,wilcox_test_func(NO.IFRs.SEPV.gnomAD.v3.Mitochondrial,IFRs.SEPV.gnomAD.v3.autosomal))
results=rbind(results,wilcox_test_func(NO.IFRs.SEPV.gnomAD.v3.X,IFRs.SEPV.gnomAD.v3.autosomal))
results=rbind(results,wilcox_test_func(NO.IFRs.SEPV.gnomAD.v3.autosomal,IFRs.SEPV.gnomAD.v3.Mitochondrial))
results=rbind(results,wilcox_test_func(NO.IFRs.SEPV.gnomAD.v3.Mitochondrial,IFRs.SEPV.gnomAD.v3.Mitochondrial))
results=rbind(results,wilcox_test_func(NO.IFRs.SEPV.gnomAD.v3.X,IFRs.SEPV.gnomAD.v3.Mitochondrial))
results=rbind(results,wilcox_test_func(NO.IFRs.SEPV.gnomAD.v3.autosomal,IFRs.SEPV.gnomAD.v3.X))
results=rbind(results,wilcox_test_func(NO.IFRs.SEPV.gnomAD.v3.Mitochondrial,IFRs.SEPV.gnomAD.v3.X))
results=rbind(results,wilcox_test_func(NO.IFRs.SEPV.gnomAD.v3.X,IFRs.SEPV.gnomAD.v3.X))

# Write tsv table with the results
write.table(results,'./Output/IFRs_SEPV_analysis/SEPV.gnomAD.v3.NO.IFRs.Vs.IFRs.by.subunit.type.tsv',sep = '\t',quote = F,row.names = F)



# Plot IFR/Non-IFR x Changing/Non Changing positions distribution

# Define function to make contingency tables to be represented in mosaic plots
NCP_IFR_NIFR=function(SEPV.IFR,SEPV.No.IFR){
  df=data.frame(IFR=c(length(SEPV.IFR[SEPV.IFR==0]),length(SEPV.IFR[SEPV.IFR>0])),Non_IFR=c(length(SEPV.No.IFR[SEPV.No.IFR==0]),length(SEPV.No.IFR[SEPV.No.IFR>0])))
  colnames(df) <- c("IFRs", "Non-IFRs")
  row.names (df) = c("Non Changing", "Changing")
  return(df)
}
# Define function to give format to p.value of Fisher test in the plots 
format.p=function(p){
  if (p<0.001)return('p < 0.001***')
  if (p<0.01 & p>0.001)return('p < 0.01**')
  if (p<0.05 & p>0.01)return('p < 0.05*')
  else return(paste0("p = ",round(p, 4)))
}
pdf("./Output/IFRs_SEPV_analysis/Figure_S3.pdf", width = 12, height = 10)
layout(matrix(c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6), nr=2, byrow=T))
par(cex=0.75)
p=fisher.test(NCP_IFR_NIFR(IFRs.SEPV.msa.interspecies.autosomal,NO.IFRs.SEPV.msa.interspecies.autosomal))$p.value
mosaicplot(t(NCP_IFR_NIFR(IFRs.SEPV.msa.interspecies.autosomal,NO.IFRs.SEPV.msa.interspecies.autosomal)), 
           main = paste0("Autosomal subunits (Inter-species)","\n"," Fisher test: ",format.p(p)),
           color = TRUE, type = 'pearson')
p=fisher.test(NCP_IFR_NIFR(IFRs.SEPV.msa.interspecies.Mitochondrial,NO.IFRs.SEPV.msa.interspecies.Mitochondrial))$p.value
mosaicplot(t(NCP_IFR_NIFR(IFRs.SEPV.msa.interspecies.Mitochondrial,NO.IFRs.SEPV.msa.interspecies.Mitochondrial)),
           main = paste0("mtDNA subunits (Inter-species)","\n"," Fisher test: ",format.p(p)),
           color = TRUE)

p=fisher.test(NCP_IFR_NIFR(IFRs.SEPV.msa.interspecies.X,NO.IFRs.SEPV.msa.interspecies.X))$p.value
mosaicplot(t(NCP_IFR_NIFR(IFRs.SEPV.msa.interspecies.X,NO.IFRs.SEPV.msa.interspecies.X)),
           main = paste0("X-linked subunits (Inter-species)","\n"," Fisher test: ", format.p(p)),
           color = TRUE)

p=fisher.test(NCP_IFR_NIFR(IFRs.SEPV.gnomAD.v3.autosomal,NO.IFRs.SEPV.gnomAD.v3.autosomal))$p.value
mosaicplot(t(NCP_IFR_NIFR(IFRs.SEPV.gnomAD.v3.autosomal,NO.IFRs.SEPV.gnomAD.v3.autosomal)),
           main = paste0("Autosomal encoded subunits (gnomAD data)","\n"," Fisher test: ",format.p(p)), 
           color = TRUE)
p=fisher.test(NCP_IFR_NIFR(IFRs.SEPV.gnomAD.v3.Mitochondrial,NO.IFRs.SEPV.gnomAD.v3.Mitochondrial))$p.value
mosaicplot(t(NCP_IFR_NIFR(IFRs.SEPV.gnomAD.v3.Mitochondrial,NO.IFRs.SEPV.gnomAD.v3.Mitochondrial)),
           main = paste0("mtDNA encoded subunits (gnomAD data)","\n"," Fisher test: ",format.p(p)),
           color = TRUE)
p=fisher.test(NCP_IFR_NIFR(IFRs.SEPV.gnomAD.v3.X,NO.IFRs.SEPV.gnomAD.v3.X))$p.value
mosaicplot(t(NCP_IFR_NIFR(IFRs.SEPV.gnomAD.v3.X,NO.IFRs.SEPV.gnomAD.v3.X)),
           main = paste0("X-linked subunits (gnomAD data)","\n"," Fisher test: ",format.p(p)), 
           color = TRUE)
dev.off()


################################################
### Comparison for SEPV from MSA within IFRs ###
################################################

# Define values of SEPV.msa.interspecies from interface residues between chains, according to the type of inheritance of involved subunits

##########
# IF A-A #
##########

ind=intersect(which(p1@genome.query=='Autosomal'),which(p1@genome.target=='Autosomal'))
IFR.A.A.SEPV.msa.interspecies=unlist(p1@SEPV.msa.interspecies[ind])
names(IFR.A.A.SEPV.msa.interspecies)=paste0(p1@subunit.query[ind],':',unlist(p1@Interface.residues.query[ind]))
IFR.A.A.SEPV.msa.interspecies=IFR.A.A.SEPV.msa.interspecies[names(IFR.A.A.SEPV.msa.interspecies)]

##########
# IF M-M #
##########

ind=intersect(which(p1@genome.query=='Mitochondrial'),which(p1@genome.target=='Mitochondrial'))
IFR.M.M.SEPV.msa.interspecies=unlist(p1@SEPV.msa.interspecies[ind])
names(IFR.M.M.SEPV.msa.interspecies)=paste0(p1@subunit.query[ind],':',unlist(p1@Interface.residues.query[ind]))
IFR.M.M.SEPV.msa.interspecies=IFR.M.M.SEPV.msa.interspecies[names(IFR.M.M.SEPV.msa.interspecies)]

##########
# IF A-M #
##########

ind=intersect(which(p1@genome.query=='Autosomal'),which(p1@genome.target=='Mitochondrial'))
IFR.A.M.SEPV.msa.interspecies=unlist(p1@SEPV.msa.interspecies[ind])
names(IFR.A.M.SEPV.msa.interspecies)=paste0(p1@subunit.query[ind],':',unlist(p1@Interface.residues.query[ind]))
IFR.A.M.SEPV.msa.interspecies=IFR.A.M.SEPV.msa.interspecies[names(IFR.A.M.SEPV.msa.interspecies)]

##########
# IF M-A #
##########

ind=intersect(which(p1@genome.query=='Mitochondrial'),which(p1@genome.target=='Autosomal'))
IFR.M.A.SEPV.msa.interspecies=unlist(p1@SEPV.msa.interspecies[ind])
names(IFR.M.A.SEPV.msa.interspecies)=paste0(p1@subunit.query[ind],':',unlist(p1@Interface.residues.query[ind]))
IFR.M.A.SEPV.msa.interspecies=IFR.M.A.SEPV.msa.interspecies[names(IFR.M.A.SEPV.msa.interspecies)]

##########
# IF X-A #
##########

ind=intersect(which(p1@genome.query=='X'),which(p1@genome.target=='Autosomal'))
IFR.X.A.SEPV.msa.interspecies=unlist(p1@SEPV.msa.interspecies[ind])
names(IFR.X.A.SEPV.msa.interspecies)=paste0(p1@subunit.query[ind],':',unlist(p1@Interface.residues.query[ind]))
IFR.X.A.SEPV.msa.interspecies=IFR.X.A.SEPV.msa.interspecies[names(IFR.X.A.SEPV.msa.interspecies)]


##########
# IF A-X #
##########

ind=intersect(which(p1@genome.query=='Autosomal'),which(p1@genome.target=='X'))
IFR.A.X.SEPV.msa.interspecies=unlist(p1@SEPV.msa.interspecies[ind])
names(IFR.A.X.SEPV.msa.interspecies)=paste0(p1@subunit.query[ind],':',unlist(p1@Interface.residues.query[ind]))
IFR.A.X.SEPV.msa.interspecies=IFR.A.X.SEPV.msa.interspecies[names(IFR.A.X.SEPV.msa.interspecies)]

##########
# IF X-M #
##########

ind=intersect(which(p1@genome.query=='X'),which(p1@genome.target=='Mitochondrial'))
IFR.X.M.SEPV.msa.interspecies=unlist(p1@SEPV.msa.interspecies[ind])
names(IFR.X.M.SEPV.msa.interspecies)=paste0(p1@subunit.query[ind],':',unlist(p1@Interface.residues.query[ind]))
IFR.X.M.SEPV.msa.interspecies=IFR.X.M.SEPV.msa.interspecies[names(IFR.X.M.SEPV.msa.interspecies)]

##########
# IF M-X #
##########

ind=intersect(which(p1@genome.query=='Mitochondrial'),which(p1@genome.target=='X'))
IFR.M.X.SEPV.msa.interspecies=unlist(p1@SEPV.msa.interspecies[ind])
names(IFR.M.X.SEPV.msa.interspecies)=paste0(p1@subunit.query[ind],':',unlist(p1@Interface.residues.query[ind]))
IFR.M.X.SEPV.msa.interspecies=IFR.M.X.SEPV.msa.interspecies[names(IFR.M.X.SEPV.msa.interspecies)]

##########
# IF X-X #
##########

# IMPORTANT: This type of IFRs does not exist


# Create a data frame with all 28 possible comparisons by IFR type for SEPV msa interspecies
results=data.frame()
results=rbind(results,wilcox_test_func(IFR.A.A.SEPV.msa.interspecies,IFR.M.A.SEPV.msa.interspecies))
results=rbind(results,wilcox_test_func(IFR.A.A.SEPV.msa.interspecies,IFR.A.M.SEPV.msa.interspecies))
results=rbind(results,wilcox_test_func(IFR.A.A.SEPV.msa.interspecies,IFR.M.M.SEPV.msa.interspecies))
results=rbind(results,wilcox_test_func(IFR.A.A.SEPV.msa.interspecies,IFR.M.X.SEPV.msa.interspecies))
results=rbind(results,wilcox_test_func(IFR.A.A.SEPV.msa.interspecies,IFR.X.M.SEPV.msa.interspecies))
results=rbind(results,wilcox_test_func(IFR.A.A.SEPV.msa.interspecies,IFR.X.A.SEPV.msa.interspecies))
results=rbind(results,wilcox_test_func(IFR.A.A.SEPV.msa.interspecies,IFR.A.X.SEPV.msa.interspecies))
results=rbind(results,wilcox_test_func(IFR.M.A.SEPV.msa.interspecies,IFR.A.M.SEPV.msa.interspecies))
results=rbind(results,wilcox_test_func(IFR.M.A.SEPV.msa.interspecies,IFR.M.M.SEPV.msa.interspecies))
results=rbind(results,wilcox_test_func(IFR.M.A.SEPV.msa.interspecies,IFR.M.X.SEPV.msa.interspecies))
results=rbind(results,wilcox_test_func(IFR.M.A.SEPV.msa.interspecies,IFR.X.M.SEPV.msa.interspecies))
results=rbind(results,wilcox_test_func(IFR.M.A.SEPV.msa.interspecies,IFR.X.A.SEPV.msa.interspecies))
results=rbind(results,wilcox_test_func(IFR.M.A.SEPV.msa.interspecies,IFR.A.X.SEPV.msa.interspecies))
results=rbind(results,wilcox_test_func(IFR.A.M.SEPV.msa.interspecies,IFR.M.M.SEPV.msa.interspecies))
results=rbind(results,wilcox_test_func(IFR.A.M.SEPV.msa.interspecies,IFR.M.X.SEPV.msa.interspecies))
results=rbind(results,wilcox_test_func(IFR.A.M.SEPV.msa.interspecies,IFR.X.M.SEPV.msa.interspecies))
results=rbind(results,wilcox_test_func(IFR.A.M.SEPV.msa.interspecies,IFR.X.A.SEPV.msa.interspecies))
results=rbind(results,wilcox_test_func(IFR.A.M.SEPV.msa.interspecies,IFR.A.X.SEPV.msa.interspecies))
results=rbind(results,wilcox_test_func(IFR.M.M.SEPV.msa.interspecies,IFR.M.X.SEPV.msa.interspecies))
results=rbind(results,wilcox_test_func(IFR.M.M.SEPV.msa.interspecies,IFR.X.M.SEPV.msa.interspecies))
results=rbind(results,wilcox_test_func(IFR.M.M.SEPV.msa.interspecies,IFR.X.A.SEPV.msa.interspecies))
results=rbind(results,wilcox_test_func(IFR.M.M.SEPV.msa.interspecies,IFR.A.X.SEPV.msa.interspecies))
results=rbind(results,wilcox_test_func(IFR.M.X.SEPV.msa.interspecies,IFR.X.M.SEPV.msa.interspecies))
results=rbind(results,wilcox_test_func(IFR.M.X.SEPV.msa.interspecies,IFR.X.A.SEPV.msa.interspecies))
results=rbind(results,wilcox_test_func(IFR.M.X.SEPV.msa.interspecies,IFR.A.X.SEPV.msa.interspecies))
results=rbind(results,wilcox_test_func(IFR.X.M.SEPV.msa.interspecies,IFR.X.A.SEPV.msa.interspecies))
results=rbind(results,wilcox_test_func(IFR.X.M.SEPV.msa.interspecies,IFR.A.X.SEPV.msa.interspecies))
results=rbind(results,wilcox_test_func(IFR.X.A.SEPV.msa.interspecies,IFR.A.X.SEPV.msa.interspecies))

# Write tsv table with the results
write.table(results,'./Output/IFRs_SEPV_analysis/SEPV.msa.interspecies.comparison.by.IFR.type.tsv',sep = '\t',quote = F,row.names = F)


######################################################
### Comparison for SEPV from gnomAD v3 within IFRs ###
######################################################

# Define values of SEPV.gnomAD.v3 from interface residues between chains, according to the type of inheritance of involved subunits

##########
# IF A-A #
##########

ind=intersect(which(p1@genome.query=='Autosomal'),which(p1@genome.target=='Autosomal'))
IFR.A.A.SEPV.gnomAD.v3=unlist(p1@SEPV.gnomAD.v3[ind])
names(IFR.A.A.SEPV.gnomAD.v3)=paste0(p1@subunit.query[ind],':',unlist(p1@Interface.residues.query[ind]))
IFR.A.A.SEPV.gnomAD.v3=IFR.A.A.SEPV.gnomAD.v3[names(IFR.A.A.SEPV.gnomAD.v3)]

##########
# IF M-M #
##########

ind=intersect(which(p1@genome.query=='Mitochondrial'),which(p1@genome.target=='Mitochondrial'))
IFR.M.M.SEPV.gnomAD.v3=unlist(p1@SEPV.gnomAD.v3[ind])
names(IFR.M.M.SEPV.gnomAD.v3)=paste0(p1@subunit.query[ind],':',unlist(p1@Interface.residues.query[ind]))
IFR.M.M.SEPV.gnomAD.v3=IFR.M.M.SEPV.gnomAD.v3[names(IFR.M.M.SEPV.gnomAD.v3)]

##########
# IF A-M #
##########

ind=intersect(which(p1@genome.query=='Autosomal'),which(p1@genome.target=='Mitochondrial'))
IFR.A.M.SEPV.gnomAD.v3=unlist(p1@SEPV.gnomAD.v3[ind])
names(IFR.A.M.SEPV.gnomAD.v3)=paste0(p1@subunit.query[ind],':',unlist(p1@Interface.residues.query[ind]))
IFR.A.M.SEPV.gnomAD.v3=IFR.A.M.SEPV.gnomAD.v3[names(IFR.A.M.SEPV.gnomAD.v3)]

##########
# IF M-A #
##########

ind=intersect(which(p1@genome.query=='Mitochondrial'),which(p1@genome.target=='Autosomal'))
IFR.M.A.SEPV.gnomAD.v3=unlist(p1@SEPV.gnomAD.v3[ind])
names(IFR.M.A.SEPV.gnomAD.v3)=paste0(p1@subunit.query[ind],':',unlist(p1@Interface.residues.query[ind]))
IFR.M.A.SEPV.gnomAD.v3=IFR.M.A.SEPV.gnomAD.v3[names(IFR.M.A.SEPV.gnomAD.v3)]

##########
# IF X-A #
##########

ind=intersect(which(p1@genome.query=='X'),which(p1@genome.target=='Autosomal'))
IFR.X.A.SEPV.gnomAD.v3=unlist(p1@SEPV.gnomAD.v3[ind])
names(IFR.X.A.SEPV.gnomAD.v3)=paste0(p1@subunit.query[ind],':',unlist(p1@Interface.residues.query[ind]))
IFR.X.A.SEPV.gnomAD.v3=IFR.X.A.SEPV.gnomAD.v3[names(IFR.X.A.SEPV.gnomAD.v3)]


##########
# IF A-X #
##########

ind=intersect(which(p1@genome.query=='Autosomal'),which(p1@genome.target=='X'))
IFR.A.X.SEPV.gnomAD.v3=unlist(p1@SEPV.gnomAD.v3[ind])
names(IFR.A.X.SEPV.gnomAD.v3)=paste0(p1@subunit.query[ind],':',unlist(p1@Interface.residues.query[ind]))
IFR.A.X.SEPV.gnomAD.v3=IFR.A.X.SEPV.gnomAD.v3[names(IFR.A.X.SEPV.gnomAD.v3)]

##########
# IF X-M #
##########

ind=intersect(which(p1@genome.query=='X'),which(p1@genome.target=='Mitochondrial'))
IFR.X.M.SEPV.gnomAD.v3=unlist(p1@SEPV.gnomAD.v3[ind])
names(IFR.X.M.SEPV.gnomAD.v3)=paste0(p1@subunit.query[ind],':',unlist(p1@Interface.residues.query[ind]))
IFR.X.M.SEPV.gnomAD.v3=IFR.X.M.SEPV.gnomAD.v3[names(IFR.X.M.SEPV.gnomAD.v3)]

##########
# IF M-X #
##########

ind=intersect(which(p1@genome.query=='Mitochondrial'),which(p1@genome.target=='X'))
IFR.M.X.SEPV.gnomAD.v3=unlist(p1@SEPV.gnomAD.v3[ind])
names(IFR.M.X.SEPV.gnomAD.v3)=paste0(p1@subunit.query[ind],':',unlist(p1@Interface.residues.query[ind]))
IFR.M.X.SEPV.gnomAD.v3=IFR.M.X.SEPV.gnomAD.v3[names(IFR.M.X.SEPV.gnomAD.v3)]

##########
# IF X-X #
##########

# IMPORTANT: This type of IFR does not exist

# Create a data frame with all 28 possible comparisons by IFR type for SEPV gnomAD.v3 values
results=data.frame()
results=rbind(results,wilcox_test_func(IFR.A.A.SEPV.gnomAD.v3,IFR.M.A.SEPV.gnomAD.v3))
results=rbind(results,wilcox_test_func(IFR.A.A.SEPV.gnomAD.v3,IFR.A.M.SEPV.gnomAD.v3))
results=rbind(results,wilcox_test_func(IFR.A.A.SEPV.gnomAD.v3,IFR.M.M.SEPV.gnomAD.v3))
results=rbind(results,wilcox_test_func(IFR.A.A.SEPV.gnomAD.v3,IFR.M.X.SEPV.gnomAD.v3))
results=rbind(results,wilcox_test_func(IFR.A.A.SEPV.gnomAD.v3,IFR.X.M.SEPV.gnomAD.v3))
results=rbind(results,wilcox_test_func(IFR.A.A.SEPV.gnomAD.v3,IFR.X.A.SEPV.gnomAD.v3))
results=rbind(results,wilcox_test_func(IFR.A.A.SEPV.gnomAD.v3,IFR.A.X.SEPV.gnomAD.v3))
results=rbind(results,wilcox_test_func(IFR.M.A.SEPV.gnomAD.v3,IFR.A.M.SEPV.gnomAD.v3))
results=rbind(results,wilcox_test_func(IFR.M.A.SEPV.gnomAD.v3,IFR.M.M.SEPV.gnomAD.v3))
results=rbind(results,wilcox_test_func(IFR.M.A.SEPV.gnomAD.v3,IFR.M.X.SEPV.gnomAD.v3))
results=rbind(results,wilcox_test_func(IFR.M.A.SEPV.gnomAD.v3,IFR.X.M.SEPV.gnomAD.v3))
results=rbind(results,wilcox_test_func(IFR.M.A.SEPV.gnomAD.v3,IFR.X.A.SEPV.gnomAD.v3))
results=rbind(results,wilcox_test_func(IFR.M.A.SEPV.gnomAD.v3,IFR.A.X.SEPV.gnomAD.v3))
results=rbind(results,wilcox_test_func(IFR.A.M.SEPV.gnomAD.v3,IFR.M.M.SEPV.gnomAD.v3))
results=rbind(results,wilcox_test_func(IFR.A.M.SEPV.gnomAD.v3,IFR.M.X.SEPV.gnomAD.v3))
results=rbind(results,wilcox_test_func(IFR.A.M.SEPV.gnomAD.v3,IFR.X.M.SEPV.gnomAD.v3))
results=rbind(results,wilcox_test_func(IFR.A.M.SEPV.gnomAD.v3,IFR.X.A.SEPV.gnomAD.v3))
results=rbind(results,wilcox_test_func(IFR.A.M.SEPV.gnomAD.v3,IFR.A.X.SEPV.gnomAD.v3))
results=rbind(results,wilcox_test_func(IFR.M.M.SEPV.gnomAD.v3,IFR.M.X.SEPV.gnomAD.v3))
results=rbind(results,wilcox_test_func(IFR.M.M.SEPV.gnomAD.v3,IFR.X.M.SEPV.gnomAD.v3))
results=rbind(results,wilcox_test_func(IFR.M.M.SEPV.gnomAD.v3,IFR.X.A.SEPV.gnomAD.v3))
results=rbind(results,wilcox_test_func(IFR.M.M.SEPV.gnomAD.v3,IFR.A.X.SEPV.gnomAD.v3))
results=rbind(results,wilcox_test_func(IFR.M.X.SEPV.gnomAD.v3,IFR.X.M.SEPV.gnomAD.v3))
results=rbind(results,wilcox_test_func(IFR.M.X.SEPV.gnomAD.v3,IFR.X.A.SEPV.gnomAD.v3))
results=rbind(results,wilcox_test_func(IFR.M.X.SEPV.gnomAD.v3,IFR.A.X.SEPV.gnomAD.v3))
results=rbind(results,wilcox_test_func(IFR.X.M.SEPV.gnomAD.v3,IFR.X.A.SEPV.gnomAD.v3))
results=rbind(results,wilcox_test_func(IFR.X.M.SEPV.gnomAD.v3,IFR.A.X.SEPV.gnomAD.v3))
results=rbind(results,wilcox_test_func(IFR.X.A.SEPV.gnomAD.v3,IFR.A.X.SEPV.gnomAD.v3))

# Write tsv table with the results
write.table(results,'./Output/IFRs_SEPV_analysis/SEPV.gnomAD.v3.comparison.by.IFR.type.tsv',sep = '\t',quote = F,row.names = F)

