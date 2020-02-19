########################################################################
#Script for tissue Specific analysis of Tiling array Expression Data
#Code developed by Elan Ness-Cohn <elanness-cohn2017@u.northwestern.edu>
########################################################################


#=====================================================================
#Dependancy Definition
#=====================================================================

library(ggplot2)
library(plyr)
library(gplots)


#=====================================================================
#Function Definition
#=====================================================================

is.not.null <- function(x) ! is.null(dim(x))

norm_samp_to_ref <- function(x,samples, refer){
  ind_S <- (which(tub_tissue_splits[[samples]]$gene %in% x))
  ind_R <- (which(tub_tissue_splits[[refer]]$gene %in% x))
  constant <- tub_tissue_splits[[refer]]$intensity[ind_R]
  tub_tissue_splits[[samples]]$intensity[ind_S]/constant
}



#=====================================================================
#Function Definition
#=====================================================================
#If set to True will create a PDF Output of all the graphs in the results folder
makePDFout = T
#If set to True will normalize the Samples Intensity's to the Reference sample
normalizeToReference = T


#=====================================================================
#Load Data and Process Data For Plotting
#=====================================================================

#Sets the working directory to /Projects Folder
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #/Projects/scripts
setwd('..') #/Projects

#Load Expression data for all 16 worm tublin genes
tub_tissue_data <- read.csv("./data/raw/TissueSpecificExpression.csv")

#Mark the Refererence Intensities with a R(reference) and others with a S(sample)
tub_tissue_data$reference <- rep("S", nrow(tub_tissue_data))
index_of_ref <- which(grepl("reference", tub_tissue_data$description))
tub_tissue_data$reference[index_of_ref] <- "R"

#tissues types to exclude based on genotypes
exclude <- c("glp-1(q224)", "dpy-28(y1) III; him-8(e1489) IV", "wdIs47 [clh-4::3XFLAG::PAB-1 + rol-6 (su1006)]", "/+; raIs/[rol-6(SU1006)+pdpy-7::GFP]","evIs111")
index <- which(tub_tissue_data$genotype %in% exclude)
tub_tissue_data <- tub_tissue_data[-index,]

#tissue tpyes to exclude based on lack of samples
exclude <- c("L1 whole animal")
index <- which(tub_tissue_data$description %in% exclude)
tub_tissue_data <- tub_tissue_data[-index,]

#make Gene Names Upper Case
tub_tissue_data$gene <- sapply(tub_tissue_data$gene, toupper)

#create Color Pal
ancestry.colours <- c('#F2F3F4', '#222222', '#F3C300', '#875692', '#F38400', '#A1CAF1', '#BE0032', '#C2B280', '#848482', '#008856', '#E68FAC', '#0067A5', '#F99379', '#604E97', '#F6A600', '#B3446C', '#DCD300', '#882D17', '#8DB600', '#654522', '#E25822', '#2B3D26')


if(normalizeToReference){
  #split the data into a list by (Beta vs Alpha), life stage, and (Reference vs Sample)
  tub_tissue_splits<-split(tub_tissue_data, list(tub_tissue_data$group,tub_tissue_data$stage,tub_tissue_data$reference))

  #remove empty dataframe from the split list (This is as a result of excluding some tissues, but tub_tissue_data$stage is still holds all levels so empty lists are beinging added to the split)
  tub_tissue_splits <- tub_tissue_splits[sapply(tub_tissue_splits,nrow)>0]

  #Initilize pairings of Sample Intensity Reads to the Reference Genes
  num_ref <- c(1:(length(names(tub_tissue_splits))/2))
  num_samples <- c(((length(names(tub_tissue_splits))/2)+1):length(names(tub_tissue_splits)))
  ref_sample_pairs <- mapply(c,num_ref,num_samples) #matrix of ref and sample pairs

  #Normalize the Intensity Reads to the Reference Genes
  data_norm <- apply(ref_sample_pairs, 2, function(y) {
    genesToCheck <- tub_tissue_splits[[y[1]]]$gene
    sapply(genesToCheck, function(gene_i) {
      norm_samp_to_ref(gene_i, samples = y[2] , refer = y[1])
    })
  })

  #set the Sample Intesisty to the Log2 normalized Sample Intensity
  apply(ref_sample_pairs,2, function(i){
    tub_tissue_splits[[i[2]]]$intensity <<- log2(as.vector(data_norm[[i[1]]]))
  })

  #only plot the normlized values
  tub_tissue_plot <- tub_tissue_splits[num_samples]

} else {
  
  #log2 the intensity Data
  tub_tissue_data$intensity <- sapply(tub_tissue_data$intensity, log2)
  
  #split the data into a list by Beta from Alpha and by life stage
  tub_tissue_splits <-split(tub_tissue_data, list(tub_tissue_data$group,tub_tissue_data$stage))
  
  #remove empty dataframe from the split list (This is as a result of excluding some tissues, but tub_tissue_data$stage is still holds all levels so empty lists beinging added to the split)
  tub_tissue_splits <- tub_tissue_splits[sapply(tub_tissue_splits,nrow)>0]
  
  #plot all samples
  tub_tissue_plot <- tub_tissue_splits
}
  

#=====================================================================
# Plotting Data
#=====================================================================

#Create the file to output all the graphs generated in the analysis
if(makePDFout){
  pdf(file = './results/normlaized_Tublin_Tissue_Processed_Expression_Plots.pdf')
}

#Create Bar Graphs for all Tublins split by Alpha and Beta. Each one of those are then seperated by Life Stage.
lapply(tub_tissue_plot, function(tub_data)
  ggplot(data = tub_data, aes(x = description, y = intensity, fill = gene)) + geom_bar(stat = "identity", position = "dodge") +
    ggtitle("Tublin Log2 Expression Values") +
    xlab("Tissue Type") +
    ylab("Intensity Log2") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5, size = 15)) +
    scale_fill_manual(values = ancestry.colours[2:10], name = "Tublin")
)

#Create Heatmaps for all Tublins split by Alpha and Beta. Each one of those are then seperated by Life Stage.
lapply(tub_tissue_plot, function(tub_data) {
  #generate a tissue by Gene matrix
  
  tissue_matrix <- aggregate(tub_data[, c("intensity")], list(tub_data$gene, tub_data$description), mean)
  tissue_matrix_p <- daply(tissue_matrix, .(Group.1, Group.2), function(vals) vals$x)
  
  #remove any columns of NA in the matrix
  if (is.not.null(tissue_matrix_p)) {
    tissue_matrix_p <- tissue_matrix_p[, !apply(is.na(tissue_matrix_p), 2, all)]
  }
  #remove any rows of NA in the matrix
  if (is.not.null(tissue_matrix_p)) {
    tissue_matrix_p <- tissue_matrix_p[!apply(is.na(tissue_matrix_p), 1, all), ]
  }
  
  if (is.not.null(tissue_matrix_p)) {
    my_palette <-colorRampPalette(c("black", "white", "yellow"))(n = 1000)
    heatmap.2(
      tissue_matrix_p,
      density.info = "none",
      main = "Log 2 Tublin Tissue \n Expression Values",
      trace = "none",
      margins = c(15, 15),
      col = my_palette,
      Colv = F,
      dendrogram = "row",
      srtCol = 45
    )
  }
})

if(makePDFout){
  dev.off()
}