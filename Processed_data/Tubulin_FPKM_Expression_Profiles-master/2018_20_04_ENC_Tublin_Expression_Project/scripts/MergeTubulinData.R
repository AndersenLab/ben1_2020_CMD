############################################################################################
#This script will merge the Tublin Data Frames for all 16 C. Elegan Tublins into a single DF  
#with the following Features: "Name" "Project" "Life.Stage" "FPKM.value" "id" "group"
#Code developed by Elan Ness-Cohn <elanness-cohn2017@u.northwestern.edu>
############################################################################################

#Set the working directory to /Projects Folder
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..') #move from script folder back to project folder

#This contains a CSV of all the genes and locations of 
masterSheet <- read.csv("./data/raw/tubulin_filePath.csv")

#assign Data for each Gene and assign it to the [gene name]_data
for(i in 1:length(masterSheet$gene)){
  assign(paste(as.character(masterSheet$gene[i]),"data", sep= "_"), read.csv(as.character(masterSheet$data[i])))
}

#addGene name Factor to each dataset
gene.data.names <- sapply(masterSheet$gene, function(geneNames) paste(geneNames,"data", sep= "_"))

#Concatinate the Data together
gene.data.names <- sapply(masterSheet$gene, function(geneNames) paste(geneNames,"data", sep= "_"))
list.of.gene.data <- lapply(gene.data.names, function(GeneList) get(GeneList))
all.Data <- do.call("rbind", list.of.gene.data)

#add the gene name to each sample
all.Data$id <- rep(masterSheet$gene, sapply(list.of.gene.data, nrow))
all.Data$group <- rep(masterSheet$group, sapply(list.of.gene.data, nrow))
names(all.Data)
write.csv(all.Data, file = "./data/raw/tubulin_data.csv")
