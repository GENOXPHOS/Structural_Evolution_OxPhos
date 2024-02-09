# Load the required libraries

library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(optparse)
library(data.table)

# Define command line input options

opls <- list()
opls[[1]] <- make_option(c("--input","-i"), action="store", type="character", dest="input", default="None", help="File containing reads with allelic imbalance measurements and their adjusted probabilities")
opls[[2]] <- make_option(c("--exp_cond","-c"), action="store", type="character", dest="condition.file", default="None", help="File containing samples distribution by condition/s to organnize the heatmap")
opls[[3]] <- make_option(c("--bed","-b"), action="store", type="character", dest="bed", default="None", help="Bed file containing analyzed genes coordinates and gene names for annotation")
opls[[4]] <- make_option(c("--exp_name","-n"), action="store", type="character", dest="exp_description", default="None", help="Character string defining the conditions or description of the experiment")
opls[[5]] <- make_option(c("--Allele_A","-A"), action="store", type="character", dest="Allele_A", default="Allele_A", help="Name of allele A")
opls[[6]] <- make_option(c("--Allele_B","-B"), action="store", type="character", dest="Allele_B", default="Allele_B", help="Name of allele B")

opts <- OptionParser(usage = "usage: %prog [options]",
                     option_list = opls, add_help_option = TRUE,
                     prog = "Rscript_4_AIM.R", description = "", epilogue = "")

args <- parse_args(opts,
                   args = commandArgs(trailingOnly = TRUE),
                   print_help_and_exit = TRUE,
                   positional_arguments = FALSE)


# Set alleles to be analyzed

Allele_A=args$Allele_A
Allele_B=args$Allele_B


############################################################
############ AIM ANALYSIS PER GENE & PER CELL ##############
############################################################

# Define the path of dyrectory containing output files from AIM.read.finder.py

path=args$input

# Create output directory

dir.create('./AIM_output/')

# Define the list of files to be analyzed

list_files=list.files(path=path,pattern = '.output.txt')

# Create an empty data frame to be updated with AIM results per gene and per cell

All <- data.frame()

# Load bed file used to annotate reads

bed <- fread(args$bed, col.names = c("chr", "start", "end","strand", "Gene_ID", "Gene_Symbol",'Gene_type'))

# Iterates each file from input directory

for (f in list_files) {
  
  # read data frames
  
  df <- read.delim(paste0(path,f),sep = "\t", stringsAsFactors = F)
  
  # Adapt bed to fit chromosome format
  if (any(grepl('chr',df$Chrom))){
      bed$chr=paste0("chr",bed$chr)
      bed$chr[bed$chr=="chrMT"]="chrM"
  }
  
  # Annotate reads with Ensembl id and gene symbol from bed file
  
  df$Gene_ID=NA
  df$Gene_Symbol=NA
  for (i in 1:nrow(bed)){
    df$Gene_ID[df$Chrom==bed$chr[i] & df$Pos <= bed$end[i] & df$Pos >= bed$star[i]]=bed$Gene_ID[i]
    df$Gene_Symbol[df$Chrom==bed$chr[i] & df$Pos <= bed$end[i] & df$Pos >= bed$star[i]]=bed$Gene_Symbol[i]
  }
  df=df[!is.na(df$Gene_ID),]
  
  # Filters data frame based on Output errors, apply a hard filter >10 reads to analyze AIM and genotype frequency
  
  df_filtered <- df %>% group_by(Gene_ID) %>% filter(n() > 10)
  df_filtered=df_filtered %>% filter(Output != 2)

  # Add a new column with the Sample name
  
  df_filtered$Sample <- strsplit(f,split = '.',fixed = T)[[1]][1]
  genes <- unique(df_filtered$Gene_ID)
  
  # fits the model by gene and cell 
  df_genes <- do.call(rbind,sapply(genes, function(x) {
    
    # Filter the data frame according to the current gene
    
    df_gene <- df_filtered[df_filtered$Gene_ID == x, ]
    
    # Perform binomial test to assess allelic imbalance and extract the p.value
    
    v.fq=table(df_gene$Output)
    p=binom.test(x = v.fq[which.max(v.fq)], n = nrow(df_gene),p = 0.5)$p.value
    if (p>0.05) df_gene$AIM=0
    if (p<=0.05 & as.numeric(names(v.fq)[which.max(v.fq)])==0) df_gene$AIM=1
    if (p<=0.05 & as.numeric(names(v.fq)[which.max(v.fq)])==1) df_gene$AIM=-1
    df_gene$p.value=p
    
    # Returns the analyzed subset of data
    
    return(df_gene)
  }, simplify = FALSE))

  # Join the subset of data analyzed for that gene in that cell to the final data frame
  
  All <- rbind(All,df_genes)
}

# Summarize AIM data per gene per cell

summary_data=All[,c("Sample","Gene_Symbol","AIM","p.value")] %>% group_by(Sample) %>% count(Sample,Gene_Symbol,AIM,p.value)
summary_data$Adj.p=p.adjust(summary_data$p.value,method = "BH")
summary_data$AIM[summary_data$Adj.p>0.05]=0

# Write a tsv data frame with AIM results per gene per cell

write.table(summary_data[,c("Sample","Gene_Symbol","AIM","p.value","Adj.p","n")], paste0("./AIM_output/AIM.",args$exp_description,".tsv"),quote = F,row.names = F,sep = "\t")

###################################################
############ HEATMAP FOR AIM RESULTS ##############
###################################################

# Define AIM matrix to be represented in the heatmap

oxphos_gene_symbol=unique(All$Gene_Symbol[which(All$Gene_ID %in% bed$Gene_ID)])
m0=tapply(summary_data$AIM, summary_data[c("Gene_Symbol", "Sample")], mean)
m1=m0[oxphos_gene_symbol,]
colnames(m1)=colnames(m0)

# Define experimental conditions to be represented in the heatmap

if(args$condition.file=='None'){
  row_split=NULL
}
if(args$condition.file!='None'){
  ExpCond=read.delim(args$condition.file, header = T, sep = "\t",dec = ".",stringsAsFactors = F)
  ExpCond$split=sapply(c(1:nrow(ExpCond)), function(x){paste(ExpCond[x,2:ncol(ExpCond)],collapse = "+")})
  ExpCond$split=factor(ExpCond$split,levels = unique(ExpCond$split),ordered = T)
  m1=m1[,colnames(m1) %in% ExpCond$ID]
  m1=m1[,order(match(colnames(m1),ExpCond$ID))]
  row_split=ExpCond$split
}


# This piece of code is specifically adapted for OxPhos genes 

Complex2Gene=data.frame(Symbol=colnames(t(m1)))
Complex2Gene$Complex=NA
Complex2Gene$Complex[which(grepl("Nduf",Complex2Gene$Symbol,ignore.case = T)==T)]="CI"
Complex2Gene$Complex[which(grepl("mt-nd",Complex2Gene$Symbol,ignore.case = T)==T)]="CI"
Complex2Gene$Complex[which(grepl("Cox",Complex2Gene$Symbol,ignore.case = T)==T)]="CIV"
Complex2Gene$Complex[which(grepl("mt-co",Complex2Gene$Symbol,ignore.case = T)==T)]="CIV"
Complex2Gene$Complex[which(grepl("NDUFA4",Complex2Gene$Symbol,ignore.case = T)==T)]="CIV"
Complex2Gene$Complex[which(grepl("NDUFA4l2",Complex2Gene$Symbol,ignore.case = T)==T)]="CIV"
Complex2Gene$Complex[which(grepl("Sdh",Complex2Gene$Symbol,ignore.case = T)==T)]="CII"
Complex2Gene$Complex[which(grepl("Uqc",Complex2Gene$Symbol,ignore.case = T)==T)]="CIII"
Complex2Gene$Complex[which(grepl("CYC1",Complex2Gene$Symbol,ignore.case = T)==T)]="CIII"
Complex2Gene$Complex[which(grepl("atp",Complex2Gene$Symbol,ignore.case = T)==T)]="CV"
Complex2Gene$Complex[which(grepl("mt-Cytb",Complex2Gene$Symbol,ignore.case = T)==T)]="CIII"
Complex2Gene$Complex[which(grepl("Cyct",Complex2Gene$Symbol,ignore.case = T)==T)]="Cc"
Complex2Gene$Complex[which(grepl("Cycs",Complex2Gene$Symbol,ignore.case = T)==T)]="Cc"
Complex2Gene=Complex2Gene[order(match(Complex2Gene$Complex, c("Cc","CI","CII","CIII","CIV","CV"))),]

# Define colors to be assigned to biallelic or monoallelic condition

col_fun = colorRamp2(c(1, 0,-1), col=c("red","yellow","purple"))

# Define heatmap parameters specific for OxPhos genes representation

if (all(!is.na(Complex2Gene$Complex))){
  m1=m1[order(match(rownames(m1),Complex2Gene$Symbol)),]
  column_split=Complex2Gene$Complex
}

if (any(is.na(Complex2Gene$Complex))){
  column_split=NULL
}

# Represent AIM genes in Heatmap for OxPhos genes

pdf(paste0("./AIM_output/AIM.",args$exp_description,".pdf"), width = 12, height = 5)
Heatmap(t(m1), na_col = "grey",row_names_gp = gpar(fontsize = 0),
        column_names_gp = gpar(col = 'black', font = 1:3,fontsize = 7),column_title_gp =  gpar(fill = "red", col = "white", border = "red"),
        row_title_gp = gpar(fontsize = 8),row_title_rot=0, column_order =rownames(m1),row_order = ExpCond$ID, row_split = row_split,col = col_fun, column_split = column_split ,
        heatmap_legend_param = list(title = "Allele Bias", at =c(1, 0,-1),labels = c(paste0("Monoallelic ",Allele_A),"Biallelic",paste0("Monoallelic ",Allele_B)), border = "black"))
dev.off()

#############################################################
############ GENE-WISE ANALYSIS OF AIM RESULTS ##############
#############################################################


# Function to analyze gene-wise AIM trend  

GeneWiseTest=function(matrix){
  
  # Only genes with at least the 75% of cells with AIM results will be considered
  
  index=unlist(sapply(1:ncol(matrix),function(c){if(sum(table(matrix[,c]))>=nrow(matrix)*75/100){return(c)}}))
  matrix=as.matrix(matrix[,unlist(sapply(1:ncol(matrix),function(c){if(sum(table(matrix[,c]))>=nrow(matrix)*75/100){return(c)}}))])
  
  genes_cast=c()
  genes_biallelic=c()
  genes_C57=c()
  
  # Parse genes in matrix genes to determine AIM trend
  
  for (c in 1:ncol(matrix)){
    p_C57=binom.test(x = length(which(na.omit(matrix[,c])>0)),n = length(which(!is.na(matrix[,c]))),p = 0.5,alternative = 'greater')$p.value
    p_biallelic=binom.test(x = length(which(na.omit(matrix[,c])==0)),n = length(which(!is.na(matrix[,c]))),p = 0.5,alternative = 'greater')$p.value
    p_cast=binom.test(x = length(which(na.omit(matrix[,c])<0)),n = length(which(!is.na(matrix[,c]))),p = 0.5,alternative = 'greater')$p.value
    if (p_C57<=0.05) genes_C57=c(genes_C57,colnames(matrix)[c])
    if (p_biallelic<=0.05) genes_biallelic=c(genes_biallelic,colnames(matrix)[c])
    if (p_cast<=0.05) genes_cast=c(genes_cast,colnames(matrix)[c])
  } 
  
  # Define variables for output data frame
  
  AIM=c(paste0("Monoallelic ",Allele_A),'Biallelic',paste0("Monoallelic ",Allele_B))
  genes=c(ifelse(length(genes_C57)>0, yes = paste(genes_C57,collapse = ', '),no = NA),
          ifelse(length(genes_biallelic)>0, yes = paste(genes_biallelic,collapse = ', '),no = NA),
          ifelse(length(genes_cast)>0, yes = paste(genes_cast,collapse = ', '),no = NA))
  return(data.frame(AIM=AIM,Genes=genes))
}

# Gene-wise trend for AIM analysis

write.table(GeneWiseTest(t(m1)), paste0("./AIM_output/AIM.GeneWise.results.",args$exp_description,".tsv"),quote = F,row.names = F,sep = "\t")
