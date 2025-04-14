library("precrec")
library("ROCR")

# Define a function to perform min-max normalization
min_max_norm = function(x) {(x - min(x)) / (max(x) - min(x))}

# ConScore calculation from SEPV values (vertebrates, mammals, primates and humans)
df=read.delim("/path/to/input_directory/OxPhos.SEPV.values.tsv", stringsAsFactors = F)

df$WHV=rank(df$SEPV.vertebrates)/nrow(df)
df$WHM=rank(df$SEPV.mammals)/nrow(df)
df$WHP=rank(df$SEPV.primates)/nrow(df)
df$WHH=rank(df$SEPV.human)/nrow(df)
df$ConScore=1-min_max_norm(df$WHH+df$WHP+df$WHM+df$WHV)

# Load AlphaMissense hotspots values for OxPhos subunits positions
AlphaMissense=read.delim("/path/to/input_directory/OxPhos_AM_hotspot_data.tsv",stringsAsFactors = F)

# Load labeled dataset of pathogenic (1) and benign (0) OxPhos missense variants 
clinvar=read.delim("/path/to/input_directory/Clinvar_OxPhos_4_benchmarcking.tsv",sep = "\t",stringsAsFactors = F)

# Add ConScore and AlphaMissense hotspots values
clinvar$ConScore=unlist(sapply(paste0(clinvar$SYMBOL,":",clinvar$Protein_position), function(x){
  df$ConScore[paste0(df$Gene_symbol,":",df$Prot.pos)==x]
}))
clinvar$AlphaMissenseHotspot=unlist(sapply(paste0(clinvar$SYMBOL,":",clinvar$Protein_position), function(x){
  AlphaMissense$mean.AM.score[paste0(AlphaMissense$Gene.name,":",AlphaMissense$position)==x]
}))

# Calculates and print ROC and Precision-Recall curves for both scores
evalmod(scores = clinvar$ConScore, labels = clinvar$Label, mode = "rocprc")
evalmod(scores = clinvar$AlphaMissenseHotspot, labels = clinvar$Label, mode = "rocprc")

# Save Clinvar dataset with labeles OxPhos missense variantes, annotated with ConScore and 
# AlphaMissense hotspots values in a tsv file
write.table(clinvar,"/path/to/output_directory/Table_S2.tsv", sep = "\t",row.names = F, quote = F)


# Save ConScore and AlphaMissense hotspots perfomance as ROC and PR curves in pdf format
pdf("/path/to/output_directory/Figure_4AB.pdf", width = 10, height = 5)

# Define the layout for the plotting area
layout(matrix(c(1,1,1,1,0,2,2,2,2), nr=1, byrow=T))

# Set graphical parameters
par(cex=1)

# Create prediction and performance objects to calculate ROC curve for ConScore
pred.ConScore=prediction(clinvar$ConScore, clinvar$Label)
perf.ConScore=performance( pred.ConScore, "tpr", "fpr" )

# Plot ROC curve
plot(perf.ConScore, col= "black")

# Create prediction and performance objects to calculate ROC curve for AlphaMissense Hotspot
pred.AM=prediction(clinvar$AlphaMissenseHotspot, clinvar$Label)
perf.AM=performance(pred.AM, "tpr", "fpr" )

# Add ROC curve plot to ConScore plot
plot(perf.AM, add = TRUE, col= "blue")


legend(x = "bottomright",  bty = "n",cex=0.75,
       legend = c(paste0("ConScore, AUC=", round(roc(clinvar$Label, clinvar$ConScore)$auc,2)),
                  paste0("AlphaMissense hotspot, AUC=", round(roc(clinvar$Label, clinvar$AlphaMissenseHotspot)$auc,2))),
       fill = c("black","blue"))


mtext("A", adj = -0.15, cex=1.5)

# Create prediction and performance objects to calculate PR curve for ConScore
perf.ConScore=performance( pred.ConScore, "prec", "rec" )

# Plot PR curve
plot(perf.ConScore, col= "black", xlim = c(0, 1), ylim = c(0, 1))

# Create performance object to calculate PR curve for AlphaMissense Hotspot
perf.AM=performance(pred.AM, "prec", "rec" )

# Plot PR curve
plot(perf.AM, add = TRUE, col= "blue", xlim = c(0, 1), ylim = c(0, 1))

# Get area under the PR-curves values
sscurves_ConScore =evalmod(scores = clinvar$ConScore, labels = clinvar$Label, mode = "rocprc")
sscurves_AM = evalmod(scores = clinvar$AlphaMissenseHotspot, labels = clinvar$Label, mode = "rocprc")

legend(x = "bottomleft",  bty = "n",cex=0.75,
       legend = c(paste0("ConScore PR-AUC=",round(precrec::auc(sscurves_ConScore)[2,4],2)),
                  paste0("AlphaMissense hotspot PR-AUC=",round(precrec::auc(sscurves_AM)[2,4],2))),
       fill = c("black","blue"))

mtext("B", adj = -0.15, cex=1.5)
dev.off()
