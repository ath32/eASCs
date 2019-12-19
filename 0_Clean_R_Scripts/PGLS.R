library(ape)
library(geiger)
library(nlme)
library(phytools)
library(caper)

### PGLS ANALYSIS ###

setwd("/Users/atho/Documents/2_Project_FSM")
data<-read.table("9_Modeling/TXT/data.txt",header=F,sep="\t", stringsAsFactors = F)
tree<-read.nexus("9_Modeling/TREE/species_list_mlr.nwk3.tre")
names(data) <- c("Species", "ASC", "Ne", "Cellularity", 'GC', 'Gene_length', 'IGR')
rownames(data) <- data$Species
comp.data<-comparative.data(tree, data, names.col="Species", vcv.dim=2, warn.dropped=TRUE)
model<-pgls(ASC~IGR, data=comp.data, lambda = "ML")
summary(model)


