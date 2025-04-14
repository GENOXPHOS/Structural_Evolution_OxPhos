library("ggplot2")
library("ggpubr")
library("dplyr")


# Load SEPV values for each subunit position in vertebrates, mammals, primates and humans
df=read.delim("/path/to/input_directory/OxPhos.SEPV.values.tsv", stringsAsFactors = F)#df=read.delim("/path/to/input_directory/OxPhos.SEPV.values.tsv", stringsAsFactors = F)

# Adjusted the dataframe to reflect the stoichiometric composition of each respiratory complex
df2=rbind(df,
          df[df$Subunit=="ENSG00000004779",],
          df[rep(rownames(df[df$Subunit=="ENSG00000152234",]),2),],
          df[df$RC=="CIII",],
          df[rep(rownames(df[df$Subunit=="ENSG00000110955",]),2),],
          df[rep(rownames(df[df$Subunit=="ENSG00000159199",]),7),],
          df[rep(rownames(df[df$Subunit=="ENSG00000135390",]),7),],
          df[rep(rownames(df[df$Subunit=="ENSG00000154518",]),7),])


# Create a function to perform a binomial test, selecting the most appropriate alternative hypothesis
binomial_test=function(x,n,p){
  test_less=binom.test(x = x, n=n,p=p,alternative = "less")
  test_greater=binom.test(x = x, n=n,p=p,alternative = "greater")
  test_2_ties=binom.test(x = x, n=n,p=p)
  if (test_less$p.value<test_greater$p.value) results=data.frame(Hypothesis="Lower than expected",Observed=test_less$estimate,Expected=p,p.value=test_less$p.value)
  if (test_less$p.value>test_greater$p.value) results=data.frame(Hypothesis="Higher than expected",Observed=test_greater$estimate,Expected=p,p.value=test_greater$p.value)
  if (test_less$p.value>0.05 & test_greater$p.value>0.05) results=data.frame(Hypothesis="Random distribution",Observed=test_2_ties$estimate,Expected=p,p.value=test_2_ties$p.value)
  return(results)
}


# Study whether fixed positions were randomly distributed across RCs, depicted in Table S6
results=data.frame()
for (rc in unique(df2$RC)){
  results=rbind(results,data.frame(RC=rc,binomial_test(x = length(which(df2$RC==rc & df2$SEPV.human==0)), n=length(which(df2$RC==rc)),p=length(which(df2$SEPV.human==0))/nrow(df2))))
}
rownames(results)=c(1:nrow(results)) 
results$adj.p=p.adjust(results$p.value,method = "BH")
results$Hypothesis[results$adj.p>0.05]="Random distribution"

# Evaluate if fixed positions are randomly distributed among RCs
results=data.frame()
for (g in unique(df2$Gene_symbol)){
  results=rbind(results,data.frame(Gene=g,binomial_test(x = length(which(df2$Gene_symbol==g & df2$SEPV.human==0)), n=length(which(df2$Gene_symbol==g)),p=length(which(df2$SEPV.human==0))/nrow(df2))))
}
rownames(results)=c(1:nrow(results)) 
results$adj.p=p.adjust(results$p.value,method = "BH")
results$Hypothesis[results$adj.p>0.05]="Random distribution"

# Save results as Table_S3.tsv
write.table(results,"/path/to/output_directory/Table_S3.tsv",sep = "\t",row.names = F, quote = F)

# Polarplot representing gene-wise enrichment analysis in Fixed positions

# Filter only significant results
results=results[results$adj.p<0.05,]

# Annotate dataframe object and adapts it to be plotted on a polar plot
results$RC=unlist(sapply(results$Gene,function(x){return(unique(df$RC[df$Gene_symbol==x]))}))
results=results %>% mutate(angle = 90 - (as.numeric(factor(interaction(Gene,RC))) - 1) * (360 / n()))

# Custom colors
colors = c("Lower than expected" = "blue", "Higher than expected" = "red") 

# Visualize results using a polar plot
f1=ggplot(results[results$adj.p<0.05,], aes(x = interaction(Gene, RC, sep = " "), y = Observed,fill = Hypothesis)) + 
  geom_bar(stat = "identity") +  # Bars with heights determined by "Observed"
  scale_fill_manual(values = colors) +  # Apply custom colors
  coord_polar(theta = "x") +
  geom_text(data=results[results$adj.p<0.05,],aes(y = 1,angle = angle,label=Gene))+ # Apply custom angles to depict gene names
  geom_hline(yintercept =length(which(df2$SEPV.human==0))/nrow(df2), color = "black", linetype = "dashed", size = 1) + 
  labs(title = "Gene-wise enrichment analysis of Fixed Positions",
       x = "",
       y = "Observed",
       fill = "hypothesis") +  # Etiqueta de la leyenda
  theme_light() + theme(axis.text.x = element_blank()) 
f1

# Save polarplot as Figure S3
pdf("/path/to/output_directory/Figure_S3.pdf",width = 8, height = 8)
f1
dev.off()
