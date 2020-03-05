library(ape)
library(geiger)
library(nlme)
library(phytools)
library(caper)

### PGLS ANALYSIS ###

#set working directory
setwd("/Users/atho/Documents/2_Project_FSM")
#load data
data<-read.table("9_Modeling/TXT/data.txt",header=F,sep="\t", stringsAsFactors = F)
#load tree
tree<-read.nexus("9_Modeling/TREE/species_list_mlr.nwk3.tre")
#set column names of data (varies depending on the data used, of course)
names(data) <- c("Species", "ASC", "Ne", "Cellularity", 'GC', 'Gene_length', 'IGR')
#set row names to species
rownames(data) <- data$Species
#create comparative data
comp.data<-comparative.data(tree, data, names.col="Species", vcv.dim=2, warn.dropped=TRUE)
#run pgls model for the desired variables (depending on the particular test)
model<-pgls(ASC~IGR, data=comp.data, lambda = "ML")
#summarise model and assess coefficients
summary(model)


