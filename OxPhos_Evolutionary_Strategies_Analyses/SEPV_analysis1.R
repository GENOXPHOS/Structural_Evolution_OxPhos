library(ggplot2)
library(ggpubr)
library(FSA)

# Load SEPV values for each subunit position in vertebrates, mammals, primates and humans
df=read.delim('/path/to/input_directory/OxPhos.SEPV.values.tsv', stringsAsFactors = F)

# Adjusted the dataframe to reflect the stoichiometric composition of each respiratory complex
df2=rbind(df,
         df[df$Subunit=='ENSG00000004779',],
         df[rep(rownames(df[df$Subunit=='ENSG00000152234',]),2),],
         df[df$RC=='CIII',],
         df[rep(rownames(df[df$Subunit=='ENSG00000110955',]),2),],
         df[rep(rownames(df[df$Subunit=='ENSG00000159199',]),7),],
         df[rep(rownames(df[df$Subunit=='ENSG00000135390',]),7),],
         df[rep(rownames(df[df$Subunit=='ENSG00000154518',]),7),])


# Define a function that computes the relative distance as Kolmogorov-Smirnov's D
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

# Perform Kruskal-Wallis tests and Dunn's post-hoc tests with Benjamini-Hochberg correction
# for SEPV differences by RC in vertebrates, mammals, primates and human
# Print the results of Kruskal-Wallis and Dunn's tests for each group
kruskal.test(SEPV.vertebrates~RC,data = df2)
df.dunn=dunnTest(SEPV.vertebrates~RC, data = df2, method = "bh") 
df.dunn
kruskal.test(SEPV.mammals~RC,data = df2)
df.dunn=dunnTest(SEPV.mammals~RC, data = df2, method = "bh") 
df.dunn
kruskal.test(SEPV.primates~RC,data = df2)
df.dunn=dunnTest(SEPV.primates~RC, data = df2, method = "bh") 
df.dunn
kruskal.test(SEPV.human~RC,data = df2)
df.dunn=dunnTest(SEPV.human~RC, data = df2, method = "bh") 
df.dunn

set.seed(1)

# Determine 1D relative distance between RCs to be plot, in vertebrates
dfplot=Related_KS_dist(list(CI=df2$SEPV.vertebrates[df2$RC=="CI"],
                            CII=df2$SEPV.vertebrates[df2$RC=="CII"],
                            CIII=df2$SEPV.vertebrates[df2$RC=="CIII"],
                            CIV=df2$SEPV.vertebrates[df2$RC=="CIV"],
                            CV=df2$SEPV.vertebrates[df2$RC=="CV"]))

p1 = ggplot(dfplot,aes(x=KS_distance,y=0))+
  xlab("")+ggtitle("Conservation in vertebrates")+ylab("Analysis by RC")+
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


# Determine 1D relative distance between RCs to be plot, in mammals
dfplot=Related_KS_dist(list(CI=df2$SEPV.mammals[df2$RC=="CI"],
                            CII=df2$SEPV.mammals[df2$RC=="CII"],
                            CIII=df2$SEPV.mammals[df2$RC=="CIII"],
                            CIV=df2$SEPV.mammals[df2$RC=="CIV"],
                            CV=df2$SEPV.mammals[df2$RC=="CV"]))

p2 = ggplot(dfplot,aes(x=KS_distance,y=0))+xlab("")+
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


# Determine 1D relative distance between RCs to be plot, in primates
dfplot=Related_KS_dist(list(CI=df2$SEPV.primates[df2$RC=="CI"],
                            CII=df2$SEPV.primates[df2$RC=="CII"],
                            CIII=df2$SEPV.primates[df2$RC=="CIII"],
                            CIV=df2$SEPV.primates[df2$RC=="CIV"],
                            CV=df2$SEPV.primates[df2$RC=="CV"]))

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

# Determine 1D relative distance between RCs to be plot, in humans
dfplot=Related_KS_dist(list(CI=df2$SEPV.human[df2$RC=="CI"],
                            CII=df2$SEPV.human[df2$RC=="CII"],
                            CIII=df2$SEPV.human[df2$RC=="CIII"],
                            CIV=df2$SEPV.human[df2$RC=="CIV"],
                            CV=df2$SEPV.human[df2$RC=="CV"]))

p4 = ggplot(dfplot,aes(x=KS_distance,y=0))+
  xlab("")+ggtitle("Conservation in human")+
  geom_point(size = 8, fill=NA,color="pink")+
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

# Perform Kruskal-Wallis tests and Dunn's post-hoc tests with Benjamini-Hochberg correction
# for SEPV differences by Subunit types (inheritance) in vertebrates, mammals, primates and human
# Print the results of Kruskal-Wallis and Dunn's tests for each group
kruskal.test(SEPV.vertebrates~Sub.type,data = df2)
df.dunn=dunnTest(SEPV.vertebrates~Sub.type, data = df2, method = "bh") 
df.dunn
kruskal.test(SEPV.mammals~Sub.type,data = df2)
df.dunn=dunnTest(SEPV.mammals~Sub.type, data = df2, method = "bh") 
df.dunn
kruskal.test(SEPV.primates~Sub.type,data = df2)
df.dunn=dunnTest(SEPV.primates~Sub.type, data = df2, method = "bh") 
df.dunn
kruskal.test(SEPV.human~Sub.type,data = df2)
df.dunn=dunnTest(SEPV.human~Sub.type, data = df2, method = "bh") 
df.dunn

# Determine 1D relative distance between subunit types to be plot, in vertebrates
dfplot=Related_KS_dist(list(A=df2$SEPV.vertebrates[df2$Sub.type=="Autosomal"],
                            M=df2$SEPV.vertebrates[df2$Sub.type=="mtDNA encoded"],
                            X=df2$SEPV.vertebrates[df2$Sub.type=="X-linked"]))


p5 = ggplot(dfplot,aes(x=KS_distance,y=0))+
  xlab("") + ggtitle("Conservation in vertebrates") + ylab("Analysis by Suninit type") +
  geom_point(size = 8, fill=NA,color="lightblue") + 
  geom_text(aes(label = Condition), col="black")+
  theme(axis.line.y = element_blank(), axis.line.x = element_line()) +
  theme(panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  ylim(c(0,0))

p5

# Determine 1D relative distance between subunit types to be plot, in mammals
dfplot=Related_KS_dist(list(A=df2$SEPV.mammals[df2$Sub.type=="Autosomal"],
                            M=df2$SEPV.mammals[df2$Sub.type=="mtDNA encoded"],
                            X=df2$SEPV.mammals[df2$Sub.type=="X-linked"]))

p6 = ggplot(dfplot,aes(x=KS_distance,y=0))+
  xlab("") + ggtitle("Conservation in mammals") +
  geom_point(size = 8, fill=NA,color="lightgreen") + 
  geom_text(aes(label = Condition), col="black")+
  theme(axis.line.y = element_blank(), axis.line.x = element_line()) +
  theme(panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  ylim(c(0,0))

p6


# Determine 1D relative distance between subunit types to be plot, in primates
dfplot=Related_KS_dist(list(A=df2$SEPV.primates[df2$Sub.type=="Autosomal"],
                            M=df2$SEPV.primates[df2$Sub.type=="mtDNA encoded"],
                            X=df2$SEPV.primates[df2$Sub.type=="X-linked"]))

p7 = ggplot(dfplot,aes(x=-KS_distance,y=0))+
  xlab("") + ggtitle("Conservation in primates") +
  geom_point(size = 8, fill=NA,color="cyan") + 
  geom_text(aes(label = Condition), col="black")+
  theme(axis.line.y = element_blank(), axis.line.x = element_line()) +
  theme(panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  ylim(c(0,0))

p7

# Determine 1D relative distance between subunit types to be plot, in humans
dfplot=Related_KS_dist(list(A=df2$SEPV.human[df2$Sub.type=="Autosomal"],
                            M=df2$SEPV.human[df2$Sub.type=="mtDNA encoded"],
                            X=df2$SEPV.human[df2$Sub.type=="X-linked"]))


p8 = ggplot(dfplot,aes(x=-KS_distance,y=0))+
  xlab("") + ggtitle("Conservation in human") +
  geom_point(size = 8, fill=NA,color="pink") + 
  geom_text(aes(label = Condition), col="black")+
  theme(axis.line.y = element_blank(), axis.line.x = element_line()) +
  theme(panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  ylim(c(0,0))

p8

# Analyze SEPV across Subunit types (inheritance) within CI for vertebrates, mammals, primates, and humans
# Perform Kruskal-Wallis tests followed by Dunn's post-hoc tests with Benjamini-Hochberg correction
# Print the results of Kruskal-Wallis and Dunn's tests for each group
kruskal.test(SEPV.vertebrates~Sub.type,data = df2[df2$RC=="CI",])
df.dunn=dunnTest(SEPV.vertebrates~Sub.type, data = df2[df2$RC=="CI",], method = "bh") 
df.dunn
kruskal.test(SEPV.mammals~Sub.type,data = df2[df2$RC=="CI",])
df.dunn=dunnTest(SEPV.mammals~Sub.type, data = df2[df2$RC=="CI",], method = "bh") 
df.dunn
kruskal.test(SEPV.primates~Sub.type,data = df2[df2$RC=="CI",])
df.dunn=dunnTest(SEPV.primates~Sub.type, data = df2[df2$RC=="CI",], method = "bh") 
df.dunn
kruskal.test(SEPV.human~Sub.type,data = df2[df2$RC=="CI",])
df.dunn=dunnTest(SEPV.human~Sub.type, data = df2[df2$RC=="CI",], method = "bh") 
df.dunn

# Calculate 1D relative distances between subunit types within CI for plotting in vertebrates
dfplot=Related_KS_dist(list(A=df2$SEPV.vertebrates[df2$Sub.type=="Autosomal" & df2$RC == "CI"],
                            M=df2$SEPV.vertebrates[df2$Sub.type=="mtDNA encoded" & df2$RC == "CI"],
                            X=df2$SEPV.vertebrates[df2$Sub.type=="X-linked" & df2$RC == "CI"]))



p9 = ggplot(dfplot,aes(x=-KS_distance,y=0))+
  xlab("") + ggtitle("Conservation in vertebrates") + 
  ylab("Analysis by Suninit type in CI") +
  geom_point(size = 8, fill=NA,color="lightblue") + 
  geom_text(aes(label = Condition), col="black")+
  theme(axis.line.y = element_blank(), axis.line.x = element_line()) +
  theme(panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  ylim(c(0,0))

p9

# Calculate 1D relative distances between subunit types within CI for plotting in mammals
dfplot=Related_KS_dist(list(A=df2$SEPV.mammals[df2$Sub.type=="Autosomal" & df2$RC == "CI"],
                            M=df2$SEPV.mammals[df2$Sub.type=="mtDNA encoded" & df2$RC == "CI"],
                            X=df2$SEPV.mammals[df2$Sub.type=="X-linked" & df2$RC == "CI"]))

p10 = ggplot(dfplot,aes(x=-KS_distance,y=0))+
  xlab("") + ggtitle("Conservation in mammals") +
  geom_point(size = 8, fill=NA,color="lightgreen") + 
  geom_text(aes(label = Condition), col="black")+
  theme(axis.line.y = element_blank(), axis.line.x = element_line()) +
  theme(panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  ylim(c(0,0))

p10

# Calculate 1D relative distances between subunit types within CI for plotting in primates
dfplot=Related_KS_dist(list(A=df2$SEPV.primates[df2$Sub.type=="Autosomal" & df2$RC == "CI"],
                            M=df2$SEPV.primates[df2$Sub.type=="mtDNA encoded" & df2$RC == "CI"],
                            X=df2$SEPV.primates[df2$Sub.type=="X-linked" & df2$RC == "CI"]))

p11 = ggplot(dfplot,aes(x=-KS_distance,y=0))+
  xlab("") + ggtitle("Conservation in primates") +
  geom_point(size = 8, fill=NA,color="cyan") + 
  geom_text(aes(label = Condition), col="black")+
  theme(axis.line.y = element_blank(), axis.line.x = element_line()) +
  theme(panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  ylim(c(0,0))

p11

# Calculate 1D relative distances between subunit types within CI for plotting in humans
dfplot=Related_KS_dist(list(A=df2$SEPV.human[df2$Sub.type=="Autosomal" & df2$RC == "CI"],
                            M=df2$SEPV.human[df2$Sub.type=="mtDNA encoded" & df2$RC == "CI"],
                            X=df2$SEPV.human[df2$Sub.type=="X-linked" & df2$RC == "CI"]))


p12 = ggplot(dfplot,aes(x=-KS_distance,y=0))+
  xlab("") + ggtitle("Conservation in human") +
  geom_point(size = 8, fill=NA,color="pink") + 
  geom_text(aes(label = Condition), col="black")+
  theme(axis.line.y = element_blank(), axis.line.x = element_line()) +
  theme(panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  ylim(c(0,0))

p12

# Analyze SEPV across Subunit types (inheritance) within CIII2 for vertebrates, mammals, 
# primates, and humans, performing two-sample Wilcoxon tests and display results
wilcox.test(df2$SEPV.vertebrates[df2$RC=="CIII" & df2$Sub.type=="Autosomal"],df2$SEPV.vertebrates[df2$RC=="CIII" & df2$Sub.type=="mtDNA encoded"],alternative='greater')
wilcox.test(df2$SEPV.mammals[df2$RC=="CIII" & df2$Sub.type=="Autosomal"],df2$SEPV.mammals[df2$RC=="CIII" & df2$Sub.type=="mtDNA encoded"],alternative='less')
wilcox.test(df2$SEPV.primates[df2$RC=="CIII" & df2$Sub.type=="Autosomal"],df2$SEPV.primates[df2$RC=="CIII" & df2$Sub.type=="mtDNA encoded"],alternative='less')
wilcox.test(df2$SEPV.human[df2$RC=="CIII" & df2$Sub.type=="Autosomal"],df2$SEPV.human[df2$RC=="CIII" & df2$Sub.type=="mtDNA encoded"],alternative='less')

# Calculate 1D relative distances between subunit types within CIII2 for plotting in vertebrates
dfplot=Related_KS_dist(list(A=df2$SEPV.vertebrates[df2$Sub.type=="Autosomal" & df2$RC == "CIII"],
                            M=df2$SEPV.vertebrates[df2$Sub.type=="mtDNA encoded" & df2$RC == "CIII"]))



p13 = ggplot(dfplot,aes(x=KS_distance,y=0))+
  xlab("") + ggtitle("Conservation in vertebrates") + ylab("Analysis by Suninit type in CIII2") +
  geom_point(size = 8, fill=NA,color="lightblue") + geom_text(aes(label = Condition), col="black")+
  theme(axis.line.y = element_blank(), axis.line.x = element_line()) +
  theme(panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  ylim(c(0,0))

p13


# Calculate 1D relative distances between subunit types within CIII2 for plotting in mammals
dfplot=Related_KS_dist(list(A=df2$SEPV.mammals[df2$Sub.type=="Autosomal" & df2$RC == "CIII"],
                            M=df2$SEPV.mammals[df2$Sub.type=="mtDNA encoded" & df2$RC == "CIII"]))

p14 = ggplot(dfplot,aes(x=KS_distance,y=0))+
  xlab("") + ggtitle("Conservation in mammals") +
  geom_point(size = 8, fill=NA,color="lightgreen") + geom_text(aes(label = Condition), col="black")+
  theme(axis.line.y = element_blank(), axis.line.x = element_line()) +
  theme(panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  ylim(c(0,0))

p14

# Calculate 1D relative distances between subunit types within CIII2 for plotting in primates
dfplot=Related_KS_dist(list(A=df2$SEPV.primates[df2$Sub.type=="Autosomal" & df2$RC == "CIII"],
                            M=df2$SEPV.primates[df2$Sub.type=="mtDNA encoded" & df2$RC == "CIII"]))

p15 = ggplot(dfplot,aes(x=-KS_distance,y=0))+
  xlab("") + ggtitle("Conservation in primates") +
  geom_point(size = 8, fill=NA,color="cyan") + geom_text(aes(label = Condition), col="black")+
  theme(axis.line.y = element_blank(), axis.line.x = element_line()) +
  theme(panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  ylim(c(0,0))

p15


# Calculate 1D relative distances between subunit types within CIII2 for plotting in humans
dfplot=Related_KS_dist(list(A=df2$SEPV.human[df2$Sub.type=="Autosomal" & df2$RC == "CIII"],
                            M=df2$SEPV.human[df2$Sub.type=="mtDNA encoded" & df2$RC == "CIII"]))


p16 = ggplot(dfplot,aes(x=KS_distance,y=0))+
  xlab("") + ggtitle("Conservation in human") +
  geom_point(size = 8, fill=NA,color="pink") + 
  geom_text(aes(label = Condition), col="black")+
  theme(axis.line.y = element_blank(), axis.line.x = element_line()) +
  theme(panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  ylim(c(0,0))

p16

# Perform Kruskal-Wallis tests followed by Dunn's post-hoc tests with Benjamini-Hochberg correction to
# analyze SEPV across Subunit types (inheritance) within CIV for vertebrates, mammals, primates, and humans
# Print the results of Kruskal-Wallis and Dunn's tests for each group
kruskal.test(SEPV.vertebrates~Sub.type,data = df2[df2$RC=="CIV",])
df.dunn=dunnTest(SEPV.vertebrates~Sub.type, data = df2[df2$RC=="CIV",], method = "bh") 
df.dunn
kruskal.test(SEPV.mammals~Sub.type,data = df2[df2$RC=="CIV",])
df.dunn=dunnTest(SEPV.mammals~Sub.type, data = df2[df2$RC=="CIV",], method = "bh") 
df.dunn
kruskal.test(SEPV.primates~Sub.type,data = df2[df2$RC=="CIV",])
df.dunn=dunnTest(SEPV.primates~Sub.type, data = df2[df2$RC=="CIV",], method = "bh") 
df.dunn
kruskal.test(SEPV.human~Sub.type,data = df2[df2$RC=="CIV",])
df.dunn=dunnTest(SEPV.human~Sub.type, data = df2[df2$RC=="CIV",], method = "bh") 
df.dunn


# Calculate 1D relative distances between subunit types within CIV for plotting in vertebrates
dfplot=Related_KS_dist(list(A=df2$SEPV.vertebrates[df2$Sub.type=="Autosomal" & df2$RC == "CIV"],
                            M=df2$SEPV.vertebrates[df2$Sub.type=="mtDNA encoded" & df2$RC == "CIV"],
                            X=df2$SEPV.vertebrates[df2$Sub.type=="X-linked" & df2$RC == "CIV"]))



p17 = ggplot(dfplot,aes(x=KS_distance,y=0))+
  xlab("") + ggtitle("Conservation in vertebrates") + 
  ylab("Analysis by Suninit type in CIV") +
  geom_point(size = 8, fill=NA,color="lightblue") + 
  geom_text(aes(label = Condition), col="black")+
  theme(axis.line.y = element_blank(), axis.line.x = element_line()) +
  theme(panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  ylim(c(0,0))
p17


# Calculate 1D relative distances between subunit types within CIV for plotting in mammals
dfplot=Related_KS_dist(list(A=df2$SEPV.mammals[df2$Sub.type=="Autosomal" & df2$RC == "CIV"],
                            M=df2$SEPV.mammals[df2$Sub.type=="mtDNA encoded" & df2$RC == "CIV"],
                            X=df2$SEPV.mammals[df2$Sub.type=="X-linked" & df2$RC == "CIV"]))

p18 = ggplot(dfplot,aes(x=KS_distance,y=0))+
  xlab("") + ggtitle("Conservation in mammals") +
  geom_point(size = 8, fill=NA,color="lightgreen") + 
  geom_text(aes(label = Condition), col="black")+
  theme(axis.line.y = element_blank(), axis.line.x = element_line()) +
  theme(panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  ylim(c(0,0))

p18

# Calculate 1D relative distances between subunit types within CIV for plotting in primates
dfplot=Related_KS_dist(list(A=df2$SEPV.primates[df2$Sub.type=="Autosomal" & df2$RC == "CIV"],
                            M=df2$SEPV.primates[df2$Sub.type=="mtDNA encoded" & df2$RC == "CIV"],
                            X=df2$SEPV.primates[df2$Sub.type=="X-linked" & df2$RC == "CIV"]))

p19 = ggplot(dfplot,aes(x=KS_distance,y=0))+
  xlab("") + ggtitle("Conservation in primates") +
  geom_point(size = 8, fill=NA,color="cyan") + 
  geom_text(aes(label = Condition), col="black")+
  theme(axis.line.y = element_blank(), axis.line.x = element_line()) +
  theme(panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  ylim(c(0,0))

p19


# Calculate 1D relative distances between subunit types within CIV for plotting in humans
dfplot=Related_KS_dist(list(A=df2$SEPV.human[df2$Sub.type=="Autosomal" & df2$RC == "CIV"],
                            M=df2$SEPV.human[df2$Sub.type=="mtDNA encoded" & df2$RC == "CIV"],
                            X=df2$SEPV.human[df2$Sub.type=="X-linked" & df2$RC == "CIV"]))


p20 = ggplot(dfplot,aes(x=-KS_distance,y=0))+
  xlab("") + ggtitle("Conservation in human") +
  geom_point(size = 8, fill=NA,color="pink") + geom_text(aes(label = Condition), col="black")+
  theme(axis.line.y = element_blank(), axis.line.x = element_line()) +
  theme(panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  ylim(c(0,0))

p20



# Assess SEPV across Subunit types (inheritance) within CV for vertebrates, mammals, primates, and humans
# Conduct two-sample Wilcoxon tests for analysis and print outcomes
wilcox.test(df2$SEPV.vertebrates[df2$RC=="CV" & df2$Sub.type=="Autosomal"],df2$SEPV.vertebrates[df2$RC=="CV" & df2$Sub.type=="mtDNA encoded"],alternative='less')
wilcox.test(df2$SEPV.mammals[df2$RC=="CV" & df2$Sub.type=="Autosomal"],df2$SEPV.mammals[df2$RC=="CV" & df2$Sub.type=="mtDNA encoded"],alternative='less')
wilcox.test(df2$SEPV.primates[df2$RC=="CV" & df2$Sub.type=="Autosomal"],df2$SEPV.primates[df2$RC=="CV" & df2$Sub.type=="mtDNA encoded"],alternative='less')
wilcox.test(df2$SEPV.human[df2$RC=="CV" & df2$Sub.type=="Autosomal"],df2$SEPV.human[df2$RC=="CV" & df2$Sub.type=="mtDNA encoded"],alternative='less')


# Calculate 1D relative distances between subunit types within CV for plotting in vertebrates
dfplot=Related_KS_dist(list(A=df2$SEPV.vertebrates[df2$Sub.type=="Autosomal" & df2$RC == "CV"],
                            M=df2$SEPV.vertebrates[df2$Sub.type=="mtDNA encoded" & df2$RC == "CV"]))



p21 = ggplot(dfplot,aes(x=KS_distance,y=0))+
  xlab("") + ggtitle("Conservation in vertebrates") + ylab("Analysis by Suninit type in CV") +
  geom_point(size = 8, fill=NA,color="lightblue") + geom_text(aes(label = Condition), col="black")+
  theme(axis.line.y = element_blank(), axis.line.x = element_line()) +
  theme(panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  ylim(c(0,0))

p21


# Calculate 1D relative distances between subunit types within CV for plotting in mammals
dfplot=Related_KS_dist(list(A=df2$SEPV.mammals[df2$Sub.type=="Autosomal" & df2$RC == "CV"],
                            M=df2$SEPV.mammals[df2$Sub.type=="mtDNA encoded" & df2$RC == "CV"]))

p22 = ggplot(dfplot,aes(x=KS_distance,y=0))+
  xlab("") + ggtitle("Conservation in mammals") +
  geom_point(size = 8, fill=NA,color="lightgreen") + geom_text(aes(label = Condition), col="black")+
  theme(axis.line.y = element_blank(), axis.line.x = element_line()) +
  theme(panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  ylim(c(0,0))

p22

# Calculate 1D relative distances between subunit types within CV for plotting in primates
dfplot=Related_KS_dist(list(A=df2$SEPV.primates[df2$Sub.type=="Autosomal" & df2$RC == "CV"],
                            M=df2$SEPV.primates[df2$Sub.type=="mtDNA encoded" & df2$RC == "CV"]))

p23 = ggplot(dfplot,aes(x=-KS_distance,y=0))+
  xlab("") + ggtitle("Conservation in primates") +
  geom_point(size = 8, fill=NA,color="cyan") + geom_text(aes(label = Condition), col="black")+
  theme(axis.line.y = element_blank(), axis.line.x = element_line()) +
  theme(panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  ylim(c(0,0))

p23


# Calculate 1D relative distances between subunit types within CV for plotting in humans
dfplot=Related_KS_dist(list(A=df2$SEPV.human[df2$Sub.type=="Autosomal" & df2$RC == "CV"],
                            M=df2$SEPV.human[df2$Sub.type=="mtDNA encoded" & df2$RC == "CV"]))


p24 = ggplot(dfplot,aes(x=KS_distance,y=0))+
  xlab("") + ggtitle("Conservation in human") +
  geom_point(size = 8, fill=NA,color="pink") + geom_text(aes(label = Condition), col="black")+
  theme(axis.line.y = element_blank(), axis.line.x = element_line()) +
  theme(panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  ylim(c(0,0))

p24

# Save plots of 1D relative distances represented as Kolmogorov-Smirnov's D across the various assessments
pdf('/path/to/output_directory/Figure_3.pdf',height = 12,width = 8.27)
ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,ncol = 4,nrow = 6,
          labels = c("A","","","","B","","","","C","","","","D","","","","E","","","","F","","",""))
dev.off()
