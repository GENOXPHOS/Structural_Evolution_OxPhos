library(ggplot2)
library(ggpubr)


# Load SEPV values for each subunit position in vertebrates, mammals, primates and humans
df=read.delim("/path/to/input_directory/OxPhos.SEPV.values.tsv", stringsAsFactors = F)

# Adjusted the dataframe to reflect the stoichiometric composition of each respiratory complex
df2=rbind(df,
          df[df$Subunit=="ENSG00000004779",],
          df[rep(rownames(df[df$Subunit=="ENSG00000152234",]),2),],
          df[df$RC=="CIII",],
          df[rep(rownames(df[df$Subunit=="ENSG00000110955",]),2),],
          df[rep(rownames(df[df$Subunit=="ENSG00000159199",]),7),],
          df[rep(rownames(df[df$Subunit=="ENSG00000135390",]),7),],
          df[rep(rownames(df[df$Subunit=="ENSG00000154518",]),7),])


# Define a function that computes the relative distance as Kolmogorov-Smirnov"s D
# It takes as argument a list of SEPV vectors "l"
Related_KS_dist=function(l){
  
  # Get the number of vectors in the list
  n = length(l)
  
  # Initialize a matrix to store KS test statistics
  ks_matrix = matrix(0, n, n)
  
  # Compute KS test statistics for each pair of vectors
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      # Perform KS test and store the statistic in the matrix
      ks_matrix[i, j] = ks.test(l[[i]], l[[j]])$statistic
      # Symmetrically assign the statistic
      ks_matrix[j, i] = ks_matrix[i, j] 
    }
  }
  # Assign row and column names to the KS matrix
  names_vectores = names(l)
  rownames(ks_matrix) = names_vectores
  colnames(ks_matrix) = names_vectores
  
  # Perform multidimensional scaling (MDS) on the KS matrix to define relative order between vectors
  mds_1d = cmdscale(as.dist(ks_matrix), k = 1)
  mds_order = rownames(mds_1d)[order(mds_1d)]
  
  # Initialize a vector to store relative KS distances
  KS_distance=c(0)
  for (i in 2:length(l)){
    KS_distance=c(KS_distance,sum(KS_distance,ks.test(l[[mds_order[i-1]]], l[[mds_order[i]]])$statistic))
  }
  
  # Return a data frame with relative KS distances between considered elements
  return(data.frame(Condition=mds_order,KS_distance=KS_distance))
}

# Load S4 object containg information about identified IFRs
IFR.list=readRDS("/path/to/input_directory/IFRs.rds")

# Filter subunits not represented in our in silico model
df3=df2[df2$Subunit %in% IFR.list@Ens.subunit.query,]

# Define binding sites (IFRs) and non-binding sites (Non-IFRs) within the OxPhos system
df3$IFR="Non-IFR"
df3$IFR[paste0(df3$Subunit,":",df3$Prot.pos) %in% unlist(sapply(c(1:length(IFR.list@Ens.subunit.query)),function(x)paste0(IFR.list@Ens.subunit.query[x],":",IFR.list@Interface.residues.query[[x]])))]="IFR"

# Define function to performs the Wilcoxon rank-sum test (Mann-Whitney U test) on two input vectors. 
# The function tests two specific hypotheses: whether vector1 is significantly less than vector2 
# and whether vector1 is significantly greater than vector2.
wilcox_test_func = function(vector1, vector2) {
  
  # Perform Wilcoxon test for "less"
  test_less = wilcox.test(vector1, vector2, alternative = "less")
  
  # Perform Wilcoxon test for "greater"
  test_greater = wilcox.test(vector1, vector2, alternative = "greater")
  
  # Check if the p-value for "less" is significant
  if (test_less$p.value < 0.05) {
    results = data.frame(Result=paste0(deparse(substitute(vector1))," Vs ",deparse(substitute(vector2))),p_value = test_less$p.value, hypothesis = "less")
  }
  
  # Check if the p-value for "greater" is significant
  if (test_greater$p.value < 0.05) {
    results = data.frame(Result=paste0(deparse(substitute(vector1))," Vs ",deparse(substitute(vector2))), p_value = test_greater$p.value, hypothesis = "greater")
  }
  
  # If no significant hypothesis was found, return "ns"
  if (test_greater$p.value > 0.05 & test_less$p.value > 0.05) {
    results = data.frame(Result= paste0(deparse(substitute(vector1))," Vs ",deparse(substitute(vector2))),p_value = NA, hypothesis = "ns")
  }
  
  return(results)
}

# Evaluate conservation of protein binding sites compared to non-binding sites within the OxPhos system 
# in vertebrates, mammals, primates and humans
wilcox_test_func(df3$SEPV.vertebrates[df3$IFR=="IFR"],df3$SEPV.vertebrates[df3$IFR=="Non-IFR"])
wilcox_test_func(df3$SEPV.mammals[df3$IFR=="IFR"],df3$SEPV.mammals[df3$IFR=="Non-IFR"])
wilcox_test_func(df3$SEPV.primates[df3$IFR=="IFR"],df3$SEPV.primates[df3$IFR=="Non-IFR"])
wilcox_test_func(df3$SEPV.human[df3$IFR=="IFR"],df3$SEPV.human[df3$IFR=="Non-IFR"])

# Evaluate conservation of protein binding sites compared to non-binding sites only within autosomal
# OxPhos subunits in vertebrates, mammals, primates and humans
wilcox_test_func(df3$SEPV.vertebrates[df3$IFR=="IFR" & df3$Sub.type=="Autosomal"],df3$SEPV.vertebrates[df3$IFR=="Non-IFR" & df3$Sub.type=="Autosomal"])
wilcox_test_func(df3$SEPV.mammals[df3$IFR=="IFR" & df3$Sub.type=="Autosomal"],df3$SEPV.mammals[df3$IFR=="Non-IFR" & df3$Sub.type=="Autosomal"])
wilcox_test_func(df3$SEPV.primates[df3$IFR=="IFR" & df3$Sub.type=="Autosomal"],df3$SEPV.primates[df3$IFR=="Non-IFR" & df3$Sub.type=="Autosomal"])
wilcox_test_func(df3$SEPV.human[df3$IFR=="IFR" & df3$Sub.type=="Autosomal"],df3$SEPV.human[df3$IFR=="Non-IFR" & df3$Sub.type=="Autosomal"])

# Evaluate conservation of protein binding sites compared to non-binding sites only within mt-DNA encoded
# OxPhos subunits in vertebrates, mammals, primates and humans
wilcox_test_func(df3$SEPV.vertebrates[df3$IFR=="IFR" & df3$Sub.type=="mtDNA encoded"],df3$SEPV.vertebrates[df3$IFR=="Non-IFR" & df3$Sub.type=="mtDNA encoded"])
wilcox_test_func(df3$SEPV.mammals[df3$IFR=="IFR" & df3$Sub.type=="mtDNA encoded"],df3$SEPV.mammals[df3$IFR=="Non-IFR" & df3$Sub.type=="mtDNA encoded"])
wilcox_test_func(df3$SEPV.primates[df3$IFR=="IFR" & df3$Sub.type=="mtDNA encoded"],df3$SEPV.primates[df3$IFR=="Non-IFR" & df3$Sub.type=="mtDNA encoded"])
wilcox_test_func(df3$SEPV.human[df3$IFR=="IFR" & df3$Sub.type=="mtDNA encoded"],df3$SEPV.human[df3$IFR=="Non-IFR" & df3$Sub.type=="mtDNA encoded"])

# Evaluate conservation of protein binding sites compared to non-binding sites only within X-linked
# OxPhos subunits in vertebrates, mammals, primates and humans
wilcox_test_func(df3$SEPV.vertebrates[df3$IFR=="IFR" & df3$Sub.type=="X-linked"],df3$SEPV.vertebrates[df3$IFR=="Non-IFR" & df3$Sub.type=="X-linked"])
wilcox_test_func(df3$SEPV.mammals[df3$IFR=="IFR" & df3$Sub.type=="X-linked"],df3$SEPV.mammals[df3$IFR=="Non-IFR" & df3$Sub.type=="X-linked"])
wilcox_test_func(df3$SEPV.primates[df3$IFR=="IFR" & df3$Sub.type=="X-linked"],df3$SEPV.primates[df3$IFR=="Non-IFR" & df3$Sub.type=="X-linked"])
wilcox_test_func(df3$SEPV.human[df3$IFR=="IFR" & df3$Sub.type=="X-linked"],df3$SEPV.human[df3$IFR=="Non-IFR" & df3$Sub.type=="X-linked"])


# Calculate and represent 1D relative distances between IFRs and Non-IFRs for plotting in vertebrates
dfplot=Related_KS_dist(list(IFR=df3$SEPV.vertebrates[df3$IFR=="IFR"],
                            NonIFR=df3$SEPV.vertebrates[df3$IFR=="Non-IFR"]))

p1 = ggplot(dfplot,aes(x=KS_distance,y=0))+
  xlab("")+ggtitle("Conservation in vertebrates")+ylab("Analysis IFR Vs Non-IFR")+
  geom_point(size = 8,color="tomato")+
  geom_text(aes(label = Condition), col="black")+
  theme(axis.line.y = element_blank(), axis.line.x = element_line())+
  theme(panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  ylim(c(0,0))

p1


# Calculate and represent 1D relative distances between IFRs and Non-IFRs for plotting in mammals
dfplot=Related_KS_dist(list(IFR=df3$SEPV.mammals[df3$IFR=="IFR"],NonIFR=df3$SEPV.mammals[df3$IFR=="Non-IFR"]))

p2 = ggplot(dfplot,aes(x=KS_distance,y=0))+
  xlab("")+ggtitle("Conservation in mammals")+
  geom_point(size = 8,color="lightgreen")+
  geom_text(aes(label = Condition), col="black")+
  theme(axis.line.y = element_blank(), axis.line.x = element_line())+
  theme(panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  ylim(c(0,0))

p2

# Calculate and represent 1D relative distances between IFRs and Non-IFRs for plotting in primates
dfplot=Related_KS_dist(list(IFR=df3$SEPV.primates[df3$IFR=="IFR"],
                            NonIFR=df3$SEPV.primates[df3$IFR=="Non-IFR"]))

p3 = ggplot(dfplot,aes(x=KS_distance,y=0))+
  xlab("")+ggtitle("Conservation in primates")+
  geom_point(size = 8,color="cyan")+
  geom_text(aes(label = Condition), col="black")+
  theme(axis.line.y = element_blank(), axis.line.x = element_line())+
  theme(panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  ylim(c(0,0))

p3



# Calculate and represent 1D relative distances between IFRs and Non-IFRs for plotting in humans
dfplot=Related_KS_dist(list(IFR=df3$SEPV.human[df3$IFR=="IFR"],
                            NonIFR=df3$SEPV.human[df3$IFR=="Non-IFR"]))

p4 = ggplot(dfplot,aes(x=KS_distance,y=0))+
  xlab("")+ggtitle("Conservation in human")+
  geom_point(size = 8,color="purple")+
  geom_text(aes(label = Condition), col="black")+
  theme(axis.line.y = element_blank(), axis.line.x = element_line())+
  theme(panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  ylim(c(0,0))

p4



# Calculate and represent 1D relative distances between autosomal IFRs and Non-IFRs for plotting in vertebrates
dfplot=Related_KS_dist(list(IFR=df3$SEPV.vertebrates[df3$IFR=="IFR" & df3$Sub.type=="Autosomal"],
                            NonIFR=df3$SEPV.vertebrates[df3$IFR=="Non-IFR" & df3$Sub.type=="Autosomal"]))

p5 = ggplot(dfplot,aes(x=KS_distance,y=0))+
  xlab("")+ggtitle("Conservation in vertebrates")+
  ylab("Analysis IFR Vs Non-IFR\nin Autosomal subunits")+
  geom_point(size = 8,color="tomato")+
  geom_text(aes(label = Condition), col="black")+
  theme(axis.line.y = element_blank(), axis.line.x = element_line())+
  theme(panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  ylim(c(0,0))

p5


# Calculate and represent 1D relative distances between autosomal IFRs and Non-IFRs for plotting in mammals
dfplot=Related_KS_dist(list(IFR=df3$SEPV.mammals[df3$IFR=="IFR" & df3$Sub.type=="Autosomal"],
                            NonIFR=df3$SEPV.mammals[df3$IFR=="Non-IFR" & df3$Sub.type=="Autosomal"]))

p6 = ggplot(dfplot,aes(x=KS_distance,y=0))+
  xlab("")+ggtitle("Conservation in mammals")+
  geom_point(size = 8,color="lightgreen")+
  geom_text(aes(label = Condition), col="black")+
  theme(axis.line.y = element_blank(), axis.line.x = element_line())+
  theme(panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  ylim(c(0,0))

p6


# Calculate and represent 1D relative distances between autosomal IFRs and Non-IFRs for plotting in primates
dfplot=Related_KS_dist(list(IFR=df3$SEPV.primates[df3$IFR=="IFR" & df3$Sub.type=="Autosomal"],
                            NonIFR=df3$SEPV.mammals[df3$IFR=="Non-IFR" & df3$Sub.type=="Autosomal"]))

p7 = ggplot(dfplot,aes(x=KS_distance,y=0))+
  xlab("")+ggtitle("Conservation in mammals")+
  geom_point(size = 8,color="cyan")+
  geom_text(aes(label = Condition), col="black")+
  theme(axis.line.y = element_blank(), axis.line.x = element_line())+
  theme(panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  ylim(c(0,0))

p7


# Calculate and represent 1D relative distances between autosomal IFRs and Non-IFRs for plotting in humans
dfplot=Related_KS_dist(list(IFR=df3$SEPV.human[df3$IFR=="IFR" & df3$Sub.type=="Autosomal"],
                            NonIFR=df3$SEPV.human[df3$IFR=="Non-IFR" & df3$Sub.type=="Autosomal"]))

p8 = ggplot(dfplot,aes(x=KS_distance,y=0))+
  xlab("")+ggtitle("Conservation in human")+
  geom_point(size = 8,color="purple")+
  geom_text(aes(label = Condition), col="black")+
  theme(axis.line.y = element_blank(), axis.line.x = element_line())+
  theme(panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  ylim(c(0,0))

p8


# Calculate and represent 1D relative distances between mtDNA encoded IFRs and Non-IFRs for plotting in vertebrates
dfplot=Related_KS_dist(list(IFR=df3$SEPV.vertebrates[df3$IFR=="IFR" & df3$Sub.type=="mtDNA encoded"],
                            NonIFR=df3$SEPV.vertebrates[df3$IFR=="Non-IFR" & df3$Sub.type=="mtDNA encoded"]))

p9 = ggplot(dfplot,aes(x=-KS_distance,y=0))+
  xlab("")+ggtitle("Conservation in vertebrates")+ylab("Analysis IFR Vs Non-IFR\nin mtDNA encoded subunits")+
  geom_point(size = 8,color="tomato")+geom_text(aes(label = Condition), col="black")+
  theme(axis.line.y = element_blank(), axis.line.x = element_line())+
  theme(panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  ylim(c(0,0))

p9


# Calculate and represent 1D relative distances between mtDNA encoded IFRs and Non-IFRs for plotting in mammals
dfplot=Related_KS_dist(list(IFR=df3$SEPV.mammals[df3$IFR=="IFR" & df3$Sub.type=="mtDNA encoded"],
                            NonIFR=df3$SEPV.mammals[df3$IFR=="Non-IFR" & df3$Sub.type=="mtDNA encoded"]))

p10 = ggplot(dfplot,aes(x=-KS_distance,y=0))+
  xlab("")+ggtitle("Conservation in mammals")+
  geom_point(size = 8,color="lightgreen")+
  geom_text(aes(label = Condition), col="black")+
  theme(axis.line.y = element_blank(), axis.line.x = element_line())+
  theme(panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  ylim(c(0,0))

p10


# Calculate and represent 1D relative distances between mtDNA encoded IFRs and Non-IFRs for plotting in primates
dfplot=Related_KS_dist(list(IFR=df3$SEPV.primates[df3$IFR=="IFR" & df3$Sub.type=="mtDNA encoded"],
                            NonIFR=df3$SEPV.primates[df3$IFR=="Non-IFR" & df3$Sub.type=="mtDNA encoded"]))

p11 = ggplot(dfplot,aes(x=-KS_distance,y=0))+
  xlab("")+ggtitle("Conservation in primates")+
  geom_point(size = 8,color="cyan")+
  geom_text(aes(label = Condition), col="black")+
  theme(axis.line.y = element_blank(), axis.line.x = element_line())+
  theme(panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  ylim(c(0,0))

p11


# Calculate and represent 1D relative distances between mtDNA encoded IFRs and Non-IFRs for plotting in humans
dfplot=Related_KS_dist(list(IFR=df3$SEPV.human[df3$IFR=="IFR" & df3$Sub.type=="mtDNA encoded"],
                            NonIFR=df3$SEPV.human[df3$IFR=="Non-IFR" & df3$Sub.type=="mtDNA encoded"]))

p12 = ggplot(dfplot,aes(x=-KS_distance,y=0))+
  xlab("")+ggtitle("Conservation in human")+
  geom_point(size = 8,color="purple")+
  geom_text(aes(label = Condition), col="black")+
  theme(axis.line.y = element_blank(), axis.line.x = element_line())+
  theme(panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  ylim(c(0,0))

p12


# Calculate and represent 1D relative distances between X-linked IFRs and Non-IFRs for plotting in vertebrates
dfplot=Related_KS_dist(list(IFR=df3$SEPV.vertebrates[df3$IFR=="IFR" & df3$Sub.type=="X-linked"],
                            NonIFR=df3$SEPV.vertebrates[df3$IFR=="Non-IFR" & df3$Sub.type=="X-linked"]))

p13 = ggplot(dfplot,aes(x=KS_distance,y=0))+
  xlab("")+ggtitle("Conservation in vertebrates")+
  ylab("Analysis IFR Vs Non-IFR \n in X-linked subunits")+
  geom_point(size = 8,color="tomato")+
  geom_text(aes(label = Condition), col="black")+
  theme(axis.line.y = element_blank(), axis.line.x = element_line())+
  theme(panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  ylim(c(0,0))

p13


# Calculate and represent 1D relative distances between X-linked IFRs and Non-IFRs for plotting in mammals
dfplot=Related_KS_dist(list(IFR=df3$SEPV.mammals[df3$IFR=="IFR" & df3$Sub.type=="X-linked"],
                            NonIFR=df3$SEPV.mammals[df3$IFR=="Non-IFR" & df3$Sub.type=="X-linked"]))

p14 = ggplot(dfplot,aes(x=KS_distance,y=0))+
  xlab("")+ggtitle("Conservation in mammals")+
  geom_point(size = 8,color="lightgreen")+
  geom_text(aes(label = Condition), col="black")+
  theme(axis.line.y = element_blank(), axis.line.x = element_line())+
  theme(panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  ylim(c(0,0))

p14


# Calculate and represent 1D relative distances between X-linked IFRs and Non-IFRs for plotting in primates
dfplot=Related_KS_dist(list(IFR=df3$SEPV.primates[df3$IFR=="IFR" & df3$Sub.type=="X-linked"],
                            NonIFR=df3$SEPV.primates[df3$IFR=="Non-IFR" & df3$Sub.type=="X-linked"]))

p15 = ggplot(dfplot,aes(x=KS_distance,y=0))+
  xlab("")+ggtitle("Conservation in primates")+
  geom_point(size = 8,color="cyan")+
  geom_text(aes(label = Condition), col="black")+
  theme(axis.line.y = element_blank(), axis.line.x = element_line())+
  theme(panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  ylim(c(0,0))

p15


# Calculate and represent 1D relative distances between X-linked IFRs and Non-IFRs for plotting in humans
dfplot=Related_KS_dist(list(IFR=df3$SEPV.human[df3$IFR=="IFR" & df3$Sub.type=="X-linked"],
                            NonIFR=df3$SEPV.human[df3$IFR=="Non-IFR" & df3$Sub.type=="X-linked"]))

p16 = ggplot(dfplot,aes(x=KS_distance,y=0))+
  xlab("")+ggtitle("Conservation in human")+
  geom_point(size = 8,color="purple")+
  geom_text(aes(label = Condition), col="black")+
  theme(axis.line.y = element_blank(), axis.line.x = element_line())+
  theme(panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  ylim(c(0,0))

p16

# Save plots of 1D relative distances represented as Kolmogorov-Smirnov"s D across the various assessments
pdf("/path/to/output_directory/Figure_S1.pdf",height = 12,width = 8.27)
ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,ncol = 4,nrow = 4,
          labels = c("A","","","","B","","","","C","","","","D","","",""))
dev.off()

