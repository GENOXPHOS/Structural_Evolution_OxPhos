# Load libraries

library("vcfR")
library("ggplot2")

# Open vcf with missense variants annotated with VEP from the 1000 Genomes Project 
vcf = read.vcfR("/path/to/input_directory/Autosomal_RC_OxPhos_phase3_1kGP_missense.vcf.gz")

# Get genotypes
gt=vcf@gt[,-1]

# Get variant coordinates and information
variants=data.frame(vcf@fix[,c(1,2,4,5)])

# Extract information from VEP annotations
CSQ=extract.info(vcf,"CSQ")
variants$Consequence=sapply(CSQ, function(x){strsplit(x,"|",fixed=T)[[1]][2]})
variants$Gene=sapply(CSQ, function(x){strsplit(x,"|",fixed=T)[[1]][5]})
variants$SYMBOL=sapply(CSQ, function(x){strsplit(x,"|",fixed=T)[[1]][4]})
variants$Protein_position=sapply(CSQ, function(x){strsplit(x,"|",fixed=T)[[1]][15]})
variants$Amino_acids=sapply(CSQ, function(x){strsplit(x,"|",fixed=T)[[1]][16]})

# Open a table curated to focus only into those variants that are going to be part of the structure, removing those placed in the signal peptide.
df=read.delim("/path/to/input_directory/OxPhos.annotation.tsv", sep = "\t",header = T,stringsAsFactors = F)
structural.filter=sapply(1:nrow(variants), function(x){return(ifelse(variants$Protein_position[x] < df$First_AA[which(df$ensembl_gene_id == variants$Gene[x])],yes = FALSE, no = TRUE))})
variants=variants[structural.filter,]
gt=gt[structural.filter,]

# Recode genotypes
gt[gt=="0|0"]=0
gt[gt=="1|1"]=0
gt[gt=="2|2"]=0
gt[gt=="1|0"]=1
gt[gt=="0|1"]=1
gt[gt=="2|0"]=1
gt[gt=="0|2"]=1
gt[gt=="2|1"]=1
gt[gt=="1|2"]=1

# Calculate the potential number of different versions of CI
# IMPORTANT: In order to calculate all possible combinations we took into account repetition of NDUFAB1
CI=gt[variants$Gene %in% df$ensembl_gene_id[df$RC=="CI"],]
class(CI)="numeric"
CI.genes=variants$Gene[variants$Gene %in% df$ensembl_gene_id[df$RC=="CI"]]

v=unlist(sapply(unique(CI.genes), function(x){
  g.mat=t(CI[CI.genes==x,])
  if (x=="ENSG00000004779") return(rep(which(rowSums(g.mat)>0),2))
  else return(which(rowSums(g.mat)>0))
}))
v1=sapply(table(v), function(x){return(2**x)})
results.CI=table(c(rep(1,ncol(CI)-length(v1)),v1))
results.CI

# Calculate the potential number of different versions of CII
CII=gt[variants$Gene %in% df$ensembl_gene_id[df$RC=="CII"],]
class(CII)="numeric"
CII.genes=variants$Gene[variants$Gene %in% df$ensembl_gene_id[df$RC=="CII"]]
v=unlist(sapply(unique(CII.genes), function(x){return(which(rowSums(t(CII[CII.genes==x,]))>0))}))
v1=sapply(table(v), function(x){return(2**x)})
results.CII=table(c(rep(1,ncol(CII)-length(v1)),v1))
results.CII

# Calculate the potential number of different versions of CIII2
# IMPORTANT:To calculate all possible combinations in CIII we considered its dimeric nature
CIII=gt[variants$Gene %in% df$ensembl_gene_id[df$RC=="CIII"],]
class(CIII)="numeric"
CIII.genes=variants$Gene[variants$Gene %in% df$ensembl_gene_id[df$RC=="CIII"]]
v=unlist(sapply(unique(CIII.genes), function(x){return(which(rowSums(t(CIII[CIII.genes==x,]))>0))}))
v1=sapply(table(v), function(x){return(2**(2*x))})
results.CIII=table(c(rep(1,ncol(CIII)-length(v1)),v1))
results.CIII

# Calculate the potential number of different versions of CIV monomer
# IMPORTANT: In order to calculate all possible combinations we took into account the existence of paralog isoforms 
CIV=gt[variants$Gene %in% df$ensembl_gene_id[df$RC=="CIV"],]
class(CIV)="numeric"
CIV.genes=variants$Gene[variants$Gene %in% df$ensembl_gene_id[df$RC=="CIV"]]
l1=list(cox5a="ENSG00000178741",cox5b="ENSG00000135940",cox6c="ENSG00000164919",cox7c="ENSG00000127184",Sub4=c("ENSG00000131055","ENSG00000131143"),
        ndufa4=c("ENSG00000185633","ENSG00000189043"), Sub6a=c("ENSG00000111775","ENSG00000156885"),Sub6b=c("ENSG00000126267","ENSG00000160471"),
        Sub7a=c("ENSG00000112695","ENSG00000115944","ENSG00000161281"),Sub8=c("ENSG00000176340","ENSG00000187581"),Sub7b=c("ENSG00000170516","ENSG00000131174"))


CIV.combinations=do.call(expand.grid, l1)

CIV.mat=matrix(1,ncol = ncol(CIV),nrow = nrow(CIV.combinations))

for (c in 1:nrow(CIV.combinations)){
  v=unlist(sapply(CIV.combinations[c,], function(x){return(which(rowSums(t(CIV[CIV.genes==x,]))>0))}))
  for (i in unique(v)){
    CIV.mat[c,i]=2**length(v[v==i])
  }
}

results.CIV=table(colSums(CIV.mat))

# Calculate the potential number of different versions of CV monomer.
# IMPORTANT: In order to calculate all possible combinations we took into account, both the stochiometric presence of each subunit
# and the existence of alternative isoforms encoded by paralog genes.
CV=gt[variants$Gene %in% df$ensembl_gene_id[df$RC=="CV"],]
class(CV)="numeric"
CV.genes=variants$Gene[variants$Gene %in% df$ensembl_gene_id[df$RC=="CV"]]


l1=list(ATP5F1A.1="ENSG00000152234",ATP5F1A.2="ENSG00000152234",ATP5F1A.3="ENSG00000152234",ATP5F1B.1="ENSG00000110955",
        ATP5F1B.2="ENSG00000110955",ATP5F1B.3="ENSG00000110955",ATP5F1C="ENSG00000165629",ATP5F1D="ENSG00000099624",
        ATP5F1E="ENSG00000124172",ATP5PB="ENSG00000116459",ATP5IF1="ENSG00000130770", ATP5MC=c("ENSG00000135390","ENSG00000154518","ENSG00000159199"),
        ATP5PD="ENSG00000167863",ATP5ME="ENSG00000169020",ATP5PF="ENSG00000154723",ATP5MJ="ENSG00000156411",
        ATP5MK="ENSG00000173915",ATP5MF="ENSG00000241468",ATP5MG="ENSG00000167283",ATP5PO="ENSG00000241837")

CV.combinations=do.call(expand.grid, l1)
for (j in 1:7) CV.combinations=cbind(CV.combinations,CV.combinations$ATP5MC)
CV.combinations
CV.mat=matrix(1,ncol = ncol(CV),nrow = nrow(CV.combinations))

for (c in 1:nrow(CV.combinations)){
  v=unlist(sapply(CV.combinations[c,], function(x){return(which(rowSums(t(CV[CV.genes==x,]))>0))}))
  for (i in unique(v)){
    CV.mat[c,i]=2**length(v[v==i])
  }
}

results.CV=table(colSums(CV.mat))

# Plot results of individual variability by Respiratory Complex regardin 1000 GP data
versions=c(names(results.CI),names(results.CII),names(results.CIII),names(results.CIV),names(results.CV))
RC=c(rep("CI",length(results.CI)),rep("CII",length(results.CII)),rep("CIII2",length(results.CIII)),rep("CIV",length(results.CIV)),rep("CV",length(results.CV)))
N=c(as.numeric(results.CI),as.numeric(results.CII),as.numeric(results.CIII),as.numeric(results.CIV),as.numeric(results.CV))
df.plot=data.frame(RC,versions,N)
df.plot$versions=as.factor(df.plot$versions)
df.plot$versions=ordered(x = df.plot$versions,levels = unique(versions)[order(as.numeric(unique(versions)))])

g1=ggplot(df.plot,aes(x=versions,y=N,colour=RC,group=RC))+
  geom_line()+geom_point()+ylim(c(0,2504))+ylab("Ner of individuals") + xlab("Ner of alternative versions of Complexes")

# Save plot in pdf file
pdf("/path/to/output_directory/Figure_1.pdf",width = 8,height = 4)
g1
dev.off()
