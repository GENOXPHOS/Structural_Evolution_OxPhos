library(ggplot2)
library(ggpubr)
library(FSA)

# Load AlphaMissense hotspots values for OxPhos subunits positions
AlphaMissense=read.delim("/path/to/input_directory/OxPhos_AM_hotspot_data.tsv",stringsAsFactors = F)

# Define a function to perform min-max normalization
min_max_norm = function(x) {(x - min(x)) / (max(x) - min(x))}

# ConScore calculation from SEPV values (vertebrates, mammals, primates and humans)
df=read.delim("/path/to/input_directory/OxPhos.SEPV.values.tsv", stringsAsFactors = F)


df$WHV=rank(df$SEPV.vertebrates)/nrow(df)
df$WHM=rank(df$SEPV.mammals)/nrow(df)
df$WHP=rank(df$SEPV.primates)/nrow(df)
df$WHH=rank(df$SEPV.human)/nrow(df)

df$ConScore=1-min_max_norm(df$WHH+df$WHP+df$WHM+df$WHV)

# Incorporate AlphaMissense hotspot values by subunit positions
for (i in 1:nrow(df)){df$AM[i]=ifelse(paste0(df$Gene_symbol[i],':',df$Prot.pos[i]) %in% paste0(AM$Gene.name,":",AM$position),yes = AM$mean.AM.score[paste0(AM$Gene.name,":",AM$position)==paste0(df$Gene_symbol[i],':',df$Prot.pos[i])],no = NA)}

# Adjusted the dataframe to reflect the stoichiometric composition of each respiratory complex
df2=rbind(df,
          df[df$Subunit=='ENSG00000004779',],
          df[rep(rownames(df[df$Subunit=='ENSG00000152234',]),2),],
          df[df$RC=='CIII',],
          df[rep(rownames(df[df$Subunit=='ENSG00000110955',]),2),],
          df[rep(rownames(df[df$Subunit=='ENSG00000159199',]),7),],
          df[rep(rownames(df[df$Subunit=='ENSG00000135390',]),7),],
          df[rep(rownames(df[df$Subunit=='ENSG00000154518',]),7),])

# Compute and print ConScore Vs AlphaMissense hotspots correlation
cor.test(x = df2$ConScore, y = df2$AM, method = 'spearman')

# Load S4 object containg information about identified IFRs
IFR.list=readRDS("/path/to/input_directory/IFRs.rds")

# Specify IFR positions based on subunit types forming the binding site
AA=unlist(sapply(which(IFR.list@genome.query=="Autosomal" & IFR.list@genome.target=="Autosomal"),function(x)paste0(IFR.list@Ens.subunit.query[x],':',IFR.list@Interface.residues.query[[x]])))
AM=unlist(sapply(which(IFR.list@genome.query=="Autosomal" & IFR.list@genome.target=="Mitochondrial"),function(x)paste0(IFR.list@Ens.subunit.query[x],':',IFR.list@Interface.residues.query[[x]])))
MA=unlist(sapply(which(IFR.list@genome.query=="Mitochondrial" & IFR.list@genome.target=="Autosomal"),function(x)paste0(IFR.list@Ens.subunit.query[x],':',IFR.list@Interface.residues.query[[x]])))
MM=unlist(sapply(which(IFR.list@genome.query=="Mitochondrial" & IFR.list@genome.target=="Mitochondrial"),function(x)paste0(IFR.list@Ens.subunit.query[x],':',IFR.list@Interface.residues.query[[x]])))
AX=unlist(sapply(which(IFR.list@genome.query=="Autosomal" & IFR.list@genome.target=="X"),function(x)paste0(IFR.list@Ens.subunit.query[x],':',IFR.list@Interface.residues.query[[x]])))
XA=unlist(sapply(which(IFR.list@genome.query=="X" & IFR.list@genome.target=="Autosomal"),function(x)paste0(IFR.list@Ens.subunit.query[x],':',IFR.list@Interface.residues.query[[x]])))
XM=unlist(sapply(which(IFR.list@genome.query=="X" & IFR.list@genome.target=="Mitochondrial"),function(x)paste0(IFR.list@Ens.subunit.query[x],':',IFR.list@Interface.residues.query[[x]])))
MX=unlist(sapply(which(IFR.list@genome.query=="Mitochondrial" & IFR.list@genome.target=="X"),function(x)paste0(IFR.list@Ens.subunit.query[x],':',IFR.list@Interface.residues.query[[x]])))

# Create a dataframe containing SEPV, ConScore, AlphaMissense values, binding site types, 
# subunit names, and corresponding protein positions
df4=data.frame(Key=c(AA,MM,AM,MA,AX,XA,MX,XM),
               Gene_symbol=unlist(sapply(c(AA,MM,AM,MA,AX,XA,MX,XM), function(x) unique(df$Gene_symbol[paste0(df$Subunit,":",df$Prot.pos) == x]))),
               Prot.pos=unlist(sapply(c(AA,MM,AM,MA,AX,XA,MX,XM), function(x) unique(df$Prot.pos[paste0(df$Subunit,":",df$Prot.pos) == x]))),
               ConScore=unlist(sapply(c(AA,MM,AM,MA,AX,XA,MX,XM), function(x) unique(df$ConScore[paste0(df$Subunit,":",df$Prot.pos) == x]))),
               AM=unlist(sapply(c(AA,MM,AM,MA,AX,XA,MX,XM), function(x) unique(df$AM[paste0(df$Subunit,":",df$Prot.pos) == x]))),
               IFR.type=c(rep("A-A",length(AA)),rep("M-M",length(MM)),rep("A-M",length(AM)),rep("M-A",length(MA)),
                          rep("A-X",length(AX)),rep("X-A",length(XA)),rep("M-X",length(MX)),rep("X-M",length(XM))))

# Perform Kruskal-Wallis tests and Dunn's post-hoc tests with Benjamini-Hochberg correction for 
# ConScore and AlphaMissense differences by binding site types considering involved subunits
# Print the results of Dunn's tests for each group
kruskal.test(ConScore~IFR.type,data = df4)
df.dunn=dunnTest(ConScore~IFR.type, data = df4, method = "bh") 
df.dunn
kruskal.test(AM~IFR.type,data = df4)
df.dunn=dunnTest(AM~IFR.type, data = df4, method = "bh") 
df.dunn

# Create the first boxplot for ConScore differences between binding site types
g1=ggplot(df4,aes(x=IFR.type, y = ConScore, fill = IFR.type))+
  geom_boxplot()+xlab("")+
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Create the first boxplot for AlphaMissense hotspot differences between binding site types
g2=ggplot(df4,aes(x=IFR.type, y = AM, fill = IFR.type))+
  geom_boxplot()+xlab("")+ylab("AlphaMissense hotspot")+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Save both plots as a PDF file
pdf('/path/to/output_directory/Figure_4CD.pdf',width = 8, height = 5)
ggarrange(g1,g2,nrow = 1,ncol = 2, labels = c("C","D"))
dev.off()