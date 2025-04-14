library(ggplot2)
library(ggpubr)
library(FSA)
library(effectsize)


# Load AlphaMissense hotspots values for OxPhos subunits positions

PSI=read.delim("/path/to/input_directory/Table_S4.tsv",stringsAsFactors = F)

# Define a function to perform min-max normalization
min_max_norm = function(x) {(x - min(x)) / (max(x) - min(x))}

# ConScore calculation from SEPV values (vertebrates, mammals, primates and humans)
df=read.delim("/path/to/input_directory/OxPhos.SEPV.values.tsv", stringsAsFactors = F)


df$WHV=rank(df$SEPV.vertebrates)/nrow(df)
df$WHM=rank(df$SEPV.mammals)/nrow(df)
df$WHP=rank(df$SEPV.primates)/nrow(df)
df$WHH=rank(df$SEPV.human)/nrow(df)

df$ConScore=1-min_max_norm(df$WHH+df$WHP+df$WHM+df$WHV)

# Incorporate PSI values by subunit positions
for (i in 1:nrow(df)){df$PSI[i]=ifelse(paste0(df$Gene_symbol[i],':',df$Prot.pos[i]) %in% paste0(PSI$Gene,":",PSI$Res), yes = PSI$PSI[paste0(PSI$Gene,":",PSI$Res)==paste0(df$Gene_symbol[i],':',df$Prot.pos[i])],no = NA)}

# Adjusted the dataframe to reflect the stoichiometric composition of each respiratory complex
df2=rbind(df,
          df[df$Subunit=='ENSG00000004779',],
          df[rep(rownames(df[df$Subunit=='ENSG00000152234',]),2),],
          df[df$RC=='CIII',],
          df[rep(rownames(df[df$Subunit=='ENSG00000110955',]),2),],
          df[rep(rownames(df[df$Subunit=='ENSG00000159199',]),7),],
          df[rep(rownames(df[df$Subunit=='ENSG00000135390',]),7),],
          df[rep(rownames(df[df$Subunit=='ENSG00000154518',]),7),])

# Filter out subunits that are not included in our structures
df2=df2[!is.na(df2$PSI),]


# Categorized the ConScore variable for analysis
df2$Q_ConScore='Q1'
df2$Q_ConScore[df2$ConScore<=quantile(df2$ConScore,0.5) & df2$ConScore>quantile(df2$ConScore,0.25)]='Q2'
df2$Q_ConScore[df2$ConScore<quantile(df2$ConScore,0.75) & df2$ConScore>quantile(df2$ConScore,0.5)]='Q3 Non-Fixed'
df2$Q_ConScore[df2$ConScore==1]='Fixed'

df2$Q_ConScore=factor(df2$Q_ConScore,levels = c('Q1','Q2','Q3 Non-Fixed','Fixed'))

# Conduct Kruskal-Wallis tests and Dunn's post-hoc tests with Benjamini-Hochberg correction 
# to analyze PSI differences across ConScore categorized segments
# Display Kruskal-Wallis outcome and Dunn's test results for each group
kruskal.test(PSI~Q_ConScore,data = df2)
df.dunn=dunnTest(PSI~Q_ConScore,data = df2, method = "bh") 
df.dunn

# Compute effect size using Cohen's D to examine the actual differences 
# among ConScore categorized segments
d1=cohens_d(df2$PSI[df2$Q_ConScore=='Q2'], df2$PSI[df2$Q_ConScore=='Q1'])
d2=cohens_d(df2$PSI[df2$Q_ConScore=='Q3 Non-Fixed'], df2$PSI[df2$Q_ConScore=='Q2'])
d3=cohens_d(df2$PSI[df2$Q_ConScore=='Fixed'], df2$PSI[df2$Q_ConScore=='Q3 Non-Fixed'])

# Create and display a boxplot showing PSI differences among ConScore categorized segments
g1=ggplot(df2, aes(x=Q_ConScore, y=PSI, fill = Q_ConScore))+
  geom_boxplot()+ylab('Probability of Structural Impact') + xlab('') +
  theme(legend.position = 'none',axis.text.x = element_text(angle = 30, vjust = 0.75, hjust = 1))
g1

# Analyze structural resilience in relation to human constraints, categorized as Fixed (unchanging) 
# and Non-Fixed positions
df2$Human_constraints="Fixed"
df2$Human_constraints[df2$SEPV.human>0]="Non-Fixed"


# Analyze and display PSI differences between Fixed and Non-fixed positions
wilcox.test(PSI~Human_constraints,data = df2,alternative='greater')
d0=cohens_d(df2$PSI[df2$Human_constraints=='Fixed'], df2$PSI[df2$Human_constraints=='Non-Fixed'])

# Boxplot depicting Fixed Vs Non-fixed differences
g2=ggplot(df2, aes(x=Human_constraints, y=PSI, fill = Human_constraints))+
  geom_boxplot()+ylab('Probability of Structural Impact') + xlab('') +
  theme(legend.position = 'none',axis.text.x = element_text(angle = 30, vjust = 0.75, hjust = 1))
g2

# Analyze PSI across Subunit types (inheritance) 
kruskal.test(PSI~Sub.type,data = df2)
df.dunn=dunnTest(PSI~Sub.type,data = df2, method = "bh") 
df.dunn

# Compute effect size using Cohen's D
d1.2=cohens_d(df2$PSI[df2$Sub.type=='mtDNA encoded'], df2$PSI[df2$Sub.type=='Autosomal'])
d2.2=cohens_d(df2$PSI[df2$Sub.type=='mtDNA encoded'], df2$PSI[df2$Sub.type=='X-linked'])

# Boxplot illustrating variations in structural resilience across different subunit types
g3=ggplot(df2, aes(x=Sub.type, y=PSI, fill = Sub.type))+
  geom_boxplot()+ylab('Probability of Structural Impact') + xlab('') +
  theme(legend.position = 'none',axis.text.x = element_text(angle = 30, vjust = 0.75, hjust = 1))
g3

# Conduct Kruskal-Wallis tests followed by Dunn's post-hoc tests with Benjamini-Hochberg correction 
# to analyze PSI differences across RC categories and present the results
kruskal.test(PSI~RC,data = df2)
df.dunn=dunnTest(PSI~RC,data = df2, method = "bh") 
df.dunn

# Compute effect size using Cohen's D
d.CI.CII=cohens_d(df2$PSI[df2$RC=='CII'], df2$PSI[df2$RC=='CI'])
d.CI.CII
d.CI.CIII=cohens_d(df2$PSI[df2$RC=='CI'], df2$PSI[df2$RC=='CIII'])
d.CI.CIII
d.CI.CIV=cohens_d(df2$PSI[df2$RC=='CIV'], df2$PSI[df2$RC=='CI'])
d.CI.CIV
d.CI.CV=cohens_d(df2$PSI[df2$RC=='CI'], df2$PSI[df2$RC=='CV'])
d.CI.CV
d.CII.CIII=cohens_d(df2$PSI[df2$RC=='CII'], df2$PSI[df2$RC=='CIII'])
d.CII.CIII
d.CII.CIV=cohens_d(df2$PSI[df2$RC=='CII'], df2$PSI[df2$RC=='CIV'])
d.CII.CIV
d.CII.CV=cohens_d(df2$PSI[df2$RC=='CII'], df2$PSI[df2$RC=='CI'])
d.CII.CV
d.CIII.CIV=cohens_d(df2$PSI[df2$RC=='CIV'], df2$PSI[df2$RC=='CIII'])
d.CIII.CIV
d.CIII.CV=cohens_d(df2$PSI[df2$RC=='CV'], df2$PSI[df2$RC=='CIII'])
d.CIII.CV
d.CIV.CV=cohens_d(df2$PSI[df2$RC=='CIV'], df2$PSI[df2$RC=='CV'])
d.CIV.CV

# Boxplot illustrating differences in PSI between Respiratory Complexes (RCs)
g4=ggplot(df2, aes(x=RC, y=PSI, fill = RC))+
  geom_boxplot()+ylab('Probability of Structural Impact') + xlab('') +
  theme(legend.position = 'none',axis.text.x = element_text(angle = 30, vjust = 0.75, hjust = 1))
g4

# Save all boxplots as Figure_S4 in pdf format
pdf('/path/to/output_directory/Figure_S2.pdf',width = 8, height = 8)
ggarrange(g1,g2,g3,g4,ncol = 2,nrow = 2,labels = c("A","B","C","D"))
dev.off()


