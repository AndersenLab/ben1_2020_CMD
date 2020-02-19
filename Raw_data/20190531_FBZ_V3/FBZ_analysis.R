library(easysorter)
library(COPASutils)
library(dplyr)
library(tidyverse)
library(cowplot)
library(mclust)
###Read in day 1 and analyze
dir <- c("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Raw_data/20190531_FBZ_V3/")

raw <- read_data(dir)

raw_noncontam <- remove_contamination(raw)

raw_nongfp <- raw_noncontam[[1]]%>%
  mutate(row.col = paste(row,col))

mylist <-matrix(0,ncol=27,nrow = 0)

mylist <-data.frame(mylist)

colnames(mylist) <- c("date","experiment","round","assay","plate","row","col","sort","TOF","EXT","time","green","yellow","red","norm.EXT","norm.green","norm.yellow","norm.red","object","call50","stage","strain","condition","control","contamination","row.col","fit$classification")

for (x in unique(raw_nongfp$strain)){
  raw_nongfp_filt1<-raw_nongfp %>%
    filter(strain == x)
  for (con in unique(raw_nongfp_filt1$condition)){
    print(con)
    if (is.na(con)){
      next}else{
        TOF_DMSO_1 <- raw_nongfp_filt1 %>%
          filter(condition==con)%>%
          filter(EXT > 100)%>%
          filter(EXT < 2000)
        TOF_DMSO <- TOF_DMSO_1%>%
          filter(EXT < (mean(TOF_DMSO_1$EXT) + 2*sd(TOF_DMSO_1$EXT))) %>%
          filter(EXT > (mean(TOF_DMSO_1$EXT) - 2*sd(TOF_DMSO_1$EXT)))
        #fit <- Mclust(TOF_DMSO[,c(10:11,13:15)],G=2)
        #if (!is.null(fit)){
        # TOF_FULL <- cbind(TOF_DMSO, fit$classification)
        #TOF_1 <- TOF_FULL %>%
        # filter(`fit$classification`==1)
        #TOF_2 <- TOF_FULL %>%
        # filter(`fit$classification`==2)
        #if(nrow(TOF_2) > nrow(TOF_1)){
        mylist <- rbind(mylist,TOF_DMSO)
        #}else{
        # mylist <- rbind(mylist,TOF_1)
      }
  }}
N2_DMSO <- raw_nongfp %>%
  filter(condition == "DMSO")%>%
  filter(strain == "N2")
N2_DMSO_filtered <- mylist %>%
  filter(condition == "DMSO")%>%
  filter(strain == "N2")
hist(N2_DMSO$TOF,probability = TRUE,breaks = 100)
d <- density(N2_DMSO_filtered$TOF)
t <- density(N2_DMSO$TOF)
lines(t,col="green")
lines(d, col="blue")
mean_noF <- mean(N2_DMSO$TOF)
med_noF <- median(N2_DMSO$TOF)
mean_F <- mean(N2_DMSO_filtered$TOF)
med_F <- median(N2_DMSO_filtered$TOF)
abline(v=mean_noF, col="green")
abline(v=med_noF,col="green")
abline(v=med_F,col="red")
abline(v=mean_F,col="red")

raw_nongfp_1 <- mylist %>%
  filter(`fit$classification` == 2)

summedplate <- sumplate(mylist,picked_plates = TRUE, directories = FALSE, quantiles = TRUE)

summedplate_1 <- summedplate %>%
  mutate(norm.n = 100)
biopruned <- bioprune(summedplate_1)

out_pruned <- prune_outliers(biopruned)

full <- regress(out_pruned)

cols <- c("N2" = "orange", "920" = "red", "1081" = "blue", "1075" = "green", "1325" = "brown", "1327" = "yellow","N2" = "orange", "919" = "red", "1082" = "blue", "1076" = "green", "1326" = "brown", "1328" = "yellow","1139"="black","1137"="pink","1317"="gray","1097"="maroon")

out_pruned$condition <- factor(out_pruned$condition, levels = c("DMSO", "7_5", "15","30","60","120"))
x <- c("N2","920","1325","1327","1081","1076")

plot <- out_pruned %>%
  filter(trait=="q90.EXT",
         !is.na(strain),
         strain %in% x) %>%
         #condition == "30")%>%
  #strain==c("N2","882","1325","1327","1082","1097","920")) %>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "882","919","920","1081","1082","1327","1328","1325","1326","1075","1076","1137","1138","1139","1140","1317","1097","1098")))%>%
  ggplot()+
  aes(x=condition, y = phenotype)+
  #geom_jitter(aes(color =fancy_strain),width=0.1)+
  #geom_jitter(width = 0.1)+
  scale_x_discrete(labels=c("DMSO" = "DMSO", "7_5" = "7.5µM","15" = "15µM", "30" = "30µM", "60" = "60µM", "120" = "120µM"))+
  theme_bw(24)+
  scale_fill_manual(name = "Strain", labels = c("N2" = "N2","882"="Del", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V","1137"="Q131L","1317"="M257I","1097"="S145F"), values = cols)+
  theme(axis.title.x=element_blank(),
        axis.text.y=element_text(size=25),
        axis.title.y=element_text(size=30))
plot +geom_smooth(alpha=0.1,aes(group=fancy_strain, colour=fancy_strain))

