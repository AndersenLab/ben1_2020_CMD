########################################################################
#Script for Life Stage analysis of FPKM Expression Data
#Code developed by Elan Ness-Cohn <elanness-cohn2017@u.northwestern.edu>
########################################################################

#=====================================================================
#Import Libraries
#=====================================================================
library(ggplot2)
library(ggrepel)
library(gplots)
library(stringr)
library(plyr)
library(RColorBrewer)

#=====================================================================
#Function Definitions
#=====================================================================
summarySE <- function(data = NULL, measurevar, groupvars = NULL, na.rm = FALSE, conf.interval = .95, .drop = TRUE) {
    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm = FALSE) {
      if (na.rm)
        sum(!is.na(x))
      else
        length(x)
    }
    
    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(
      data,
      groupvars,
      .drop = .drop,
      .fun = function(xx, col) {
        c(
          N    = length2(xx[[col]], na.rm = na.rm),
          mean = mean   (xx[[col]], na.rm = na.rm),
          sd   = sd     (xx[[col]], na.rm = na.rm)
        )
      },
      measurevar
    )
    
    # Rename the "mean" column
    datac <- rename(datac, c("mean" = measurevar))
    
    datac$se <-
      datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
    
    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval:
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval / 2 + .5, datac$N - 1)
    datac$ci <- datac$se * ciMult
    
    return(datac)
  }

#=====================================================================
#Load Data and Process Data For Plotting
#=====================================================================

#Set the working directory to /Projects Folder
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..') #move from script folder back to project folder

#Load FPKM data for all 16 worm tublin genes
tub_FPKM_data <- read.csv("./data/raw/tubulin_data.csv")

#Load SRX file and Remove extra " and spaces, as a result of the initial download of the data
meta <- read.csv("./data/raw/control_SRX.csv")
meta$SAMPLE_ID <- str_replace_all(as.character(meta$SAMPLE_ID), '"', "")
meta$SAMPLE_ID <- str_replace_all(meta$SAMPLE_ID, ' ', "")

#log transfrom Data
tub_FPKM_data$FPKM.value <- sapply(tub_FPKM_data$FPKM.value, function(x) log2(x+1))

#replace "_" with "-" in gene name (ie. TBB_1 to TBB-1)
tub_FPKM_data$id<- sub(x = tub_FPKM_data$id, pattern = "_", replacement = "-")

#select only N2 worms
N2worms <- grep(pattern = "*.N2.*", tub_FPKM_data$Name, value = TRUE)
N2worms.index <- which(tub_FPKM_data$Name %in% N2worms)
N2_FPKM_data <- tub_FPKM_data[N2worms.index,]

#extract only experiments FPKM values whose condition are not treated with pathogen, drugs, etc. and filter data frame only to those experiments
SRX_ID <- sapply(N2_FPKM_data$Name, function(ids) sub(".*\\.", "", ids))
control_experiements <- which(SRX_ID %in% meta$SAMPLE_ID)
N2_FPKM_data_cont <- N2_FPKM_data[control_experiements,]
unique(N2_FPKM_data_cont$Life.Stage)

#get List of all Life stages in dataset
N2stages <- unique(N2_FPKM_data$Life.Stage)


#=====================================================================
#plot Minutes Post First Cleavage Data
#=====================================================================


#Create the file to output all the graphs generated in the analysis
pdf(file = './results/Tublin_FPKM_Processed_Expression_Plots.pdf')

#select only "min post first-Cleavage Ce" (mpfc) life stage
#gets names of stages at all stages between 10 minutes and 850 minutes
mpfc.Set <- grep('min post first-cleavage Ce', N2stages, value=TRUE) 
mpfc.index <- which(N2_FPKM_data$Life.Stage %in% mpfc.Set)
mpfc.N2 <- N2_FPKM_data[mpfc.index,]

#extract just the numerical time point between 10 and 850 minutes (ie. get rid of "min post first-cleavage Ce")
mpfc.N2$Life.Stage <- as.numeric(sub(pattern = "*.min post first-cleavage Ce", replacement = "", x = mpfc.N2$Life.Stage))

#aggregate the dataframe into the means of each subgrouping (time and gene)
#ie. time = 20, gene = TBA_1 -> the mean of all TBA-1 genes at timepoint 20.
mpfc.mean.matrix <- aggregate(mpfc.N2[,c("FPKM.value")], list(mpfc.N2$id,mpfc.N2$Life.Stage), mean)
mpfc.sd.matrix <- aggregate(mpfc.N2[,c("FPKM.value")], list(mpfc.N2$id,mpfc.N2$Life.Stage), sd)

#reformate the data into a matrix to be put into the plot
mpfc.matrix.plot.data <- daply(mpfc.mean.matrix, .(Group.1, Group.2), function(vals) vals$x)

#create a color palate for FPKM data
my_palette <- colorRampPalette(c("black", "white", "yellow"))(n = 1000)

#generate the heat map of the FPKM data
heatmap.2(mpfc.matrix.plot.data, 
          density.info="none",
          main = "Log 2 Tublin FPKM Expression Values in \n N2 Untreated Worms",
          xlab = "Minutes Post First Cleavage",
          trace="none",
          margins =c(10,10),
          col = my_palette,
          Colv=FALSE,
          dendrogram = "row"
          )

#plot without the common tubulins present
common_tubs <- c("TBB-1","TBB-2","TBA-1","TBA-2")
heatmap.2(mpfc.matrix.plot.data[-c(which(rownames(mpfc.matrix.plot.data)%in% common_tubs)),], 
          density.info="none",
          main = "Log 2 Tublin FPKM Expression Values in \n N2 Untreated Worms",
          xlab = "Minutes Post First Cleavage",
          trace="none",
          margins =c(10,10),
          col = my_palette,
          Colv=FALSE,
          dendrogram = "row"
)
stats.mpfc <- summarySE(mpfc.N2, measurevar="FPKM.value", groupvars=c("Life.Stage","id", "group"), na.rm = T)

#need to define. many with not enough data points
#ggplot(data=stats.mpfc, aes(x=Life.Stage, y=FPKM.value, fill = id)) + geom_line(stat = "identity", position = "dodge") + geom_errorbar(aes(ymin=FPKM.value-sd, ymax=FPKM.value+sd, na.rm = TRUE), width=.1, position=position_dodge(.1))

#Create Color Patterns
ancestry.colours <- c('#F2F3F4', '#222222', '#F3C300', '#875692', '#F38400', '#A1CAF1', '#BE0032', '#C2B280', '#848482', '#008856', '#E68FAC', '#0067A5', '#F99379', '#604E97', '#F6A600', '#B3446C', '#DCD300', '#882D17', '#8DB600', '#654522', '#E25822', '#2B3D26')

#Seperate the Alpha Tublins From the Beta Tublins
AlphaTub_FC <-stats.mpfc[which(stats.mpfc$group == "Alpha"),]
BetaTub_FC <-stats.mpfc[which(stats.mpfc$group == "Beta"),]

#plot Alpha Tublin
ggplot(data=AlphaTub_FC, aes(x=Life.Stage, y=FPKM.value, col = id)) +
  geom_line(size=1) +
  geom_errorbar(aes(ymin=FPKM.value-sd, ymax=FPKM.value+sd), width=10,size=0.5) +
  scale_colour_manual(values = ancestry.colours[2:10], name = "Alpha Tublin") +
  ylab("log2 FPKM") + 
  xlab("Minutes Post First Cleavage") +
  ggtitle("Alpha Tublin log2 FPKM Expression \n Minutes Post First Cleavage")+
  theme(plot.title = element_text(hjust = 0.5, size = 15))

#plot Beta Tublin
ggplot(data=BetaTub_FC, aes(x=Life.Stage, y=FPKM.value, col = id)) +
  geom_line(size=1) +
  geom_errorbar(aes(ymin=FPKM.value-sd, ymax=FPKM.value+sd), width=10,size=0.5) +
  scale_colour_manual(values = ancestry.colours[2:10], name = "Beta Tublin") +
  ylab("log2 FPKM") + 
  xlab("Minutes Post First Cleavage") +
  ggtitle("Beta Tublin log2 FPKM Expression \n Minutes Post First Cleavage")+
  theme(plot.title = element_text(hjust = 0.5, size = 15))


#=====================================================================
#plot the worm stages
#=====================================================================


#select only life stages
worm.stages <- c("embryo Ce" , "L1 larva Ce", "L2 larva Ce" , "L3 larva Ce" , "L4 larva Ce", "adult Ce")
worms.stages.index <- which(N2_FPKM_data$Life.Stage %in% worm.stages)
Stage.N2 <- N2_FPKM_data[worms.stages.index,]

#Convert life stages from factors to character
Stage.N2$Life.Stage <- as.character(Stage.N2$Life.Stage)

#aggregate the dataframe into the means of each subgrouping (time and gene)
#ie. time = 20, gene = TBA_1
Stage.N2.mean.matrix <- aggregate(Stage.N2[,c("FPKM.value")], list(Stage.N2$id,Stage.N2$Life.Stage), mean)
Stage.N2.sd.matrix <- aggregate(Stage.N2[,c("FPKM.value")], list(Stage.N2$id,Stage.N2$Life.Stage), sd)

#replace "_" with "-" in gene name (ie. TBB_1 to TBB-1)
Stage.N2.mean.matrix$Group.1 <- sub(x = Stage.N2.mean.matrix$Group.1, pattern = "_", replacement = "-")

#reformate the data into a matrix to be put into the plot
Stage.N2.matrix.plot.data <- daply(Stage.N2.mean.matrix, .(Group.1, Group.2), function(vals) vals$x)
Stage.N2.matrix.plot.data <-  Stage.N2.matrix.plot.data[,c(2,3,4,5,6,1)]
colnames(Stage.N2.matrix.plot.data) <- c("Embryo" , "L1", "L2" , "L3" , "L4", "Adult")

#generate the heat map of the FPKM data
heatmap.2(Stage.N2.matrix.plot.data, 
          density.info="none",
          main = "Log 2 Tublin FPKM Expression Values \n in N2 Untreated Worms",
          xlab = "Life Cycle Stage",
          trace="none",
          margins =c(10,10),
          col = my_palette,
          Colv=FALSE,
          dendrogram = "row"
)

#plot without Common Tublins
heatmap.2(Stage.N2.matrix.plot.data[-c(which(rownames(Stage.N2.matrix.plot.data) %in% common_tubs)),], 
          density.info="none",
          main = "Log 2 Tublin FPKM Expression Values in \n N2 Untreated Worms",
          xlab = "Minutes Post First Cleavage",
          trace="none",
          margins =c(10,10),
          col = my_palette,
          Colv=FALSE,
          dendrogram = "row"
)

#get Summary of the Stats for N2 worms in EE -> Adult Life stages
stats <- summarySE(Stage.N2, measurevar="FPKM.value", groupvars=c("Life.Stage","id", "group"))
stats$Life.Stage <- as.numeric(as.factor(stats$Life.Stage)) -1
stats$Life.Stage[which(stats$Life.Stage == 0)] <- 6

#Select Alpha and Beta Tublins from the list
AlphaTub_LS <-stats[which(stats$group == "Alpha"),]
BetaTub_LS <-stats[which(stats$group == "Beta"),]

#plot Alpha tublins
ggplot(data=AlphaTub_LS, aes(x=Life.Stage, y=FPKM.value, col = id)) +
  geom_line(size=1) +
  geom_errorbar(aes(ymin=FPKM.value-sd, ymax=FPKM.value+sd), width=.1,size=0.5) +
  scale_colour_manual(values = ancestry.colours[2:10], name = "Alpha Tublin") +
  ylab("log2 FPKM") + 
  xlab("Life Stage") +
  ggtitle("Alpha Tublin log2 FPKM Expression") + 
  scale_x_continuous(breaks=c(1:6), labels=c("Embryo" , "L1", "L2" , "L3" , "L4", "Adult"))+
  theme(plot.title = element_text(hjust = 0.5, size = 15))

#plot Beta tublins
ggplot(data=BetaTub_LS, aes(x=Life.Stage, y=FPKM.value, col = id)) +
  geom_line(size=1) +
  geom_errorbar(aes(ymin=FPKM.value-sd, ymax=FPKM.value+sd), width=.1,size=0.5) +
  scale_colour_manual(values = ancestry.colours[2:10], name = "Beta Tublin") +
  ylab("log2 FPKM") + 
  xlab("Life Stage") +
  ggtitle("Beta Tublin log2 FPKM Expression") + 
  scale_x_continuous(breaks=c(1:6), labels=c("Embryo" , "L1", "L2" , "L3" , "L4", "Adult"))+
  theme(plot.title = element_text(hjust = 0.5, size = 15))
  
save(BetaTub_LS, file = "~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Processed_data/Lifestagebetas.Rdata")
#turn off dev to move plots to the PDF output
dev.off()
