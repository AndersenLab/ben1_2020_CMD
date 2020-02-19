library(easysorter)
library(COPASutils)
library(dplyr)
library(tidyverse)
library(cowplot)
library(mclust)

###Read in ABZ V3 dose from 20190607
setwd("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Processed_data/")
dir <- c("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Raw_data/20190607_ABZv3_analysis/")

raw <- easysorter::read_data(dir)

raw_noncontam <- easysorter::remove_contamination(raw)


#Put into list format for Mclust analysis
raw_nongfp <- raw_noncontam[[1]]%>%
  mutate(row.col = paste(row,col))

#Cycle through each strain and condition combination and run Mclust on them using log(EXT) and log(TOF) to seperate bleach debree from worms

mylist <- NULL
for (x in unique(raw_nongfp$strain)){
  raw_nongfp_filt1<-raw_nongfp %>%
    filter(strain == x)
  for (con in unique(raw_nongfp_filt1$condition)){
    print(con)
    if (is.na(con)){
      next}else{
        TOF_DMSO_1 <- raw_nongfp_filt1 %>%
          filter(condition==con)
        TOF_DMSO <- TOF_DMSO_1%>%
          mutate(log_EXT = log(EXT))%>%
          mutate(log_TOF = log(TOF))
        fit <- Mclust(TOF_DMSO[,c(27:28)],G=2)
        if (!is.null(fit)){
          TOF_FULL <- cbind(TOF_DMSO, fit$classification)
          TOF_1 <- TOF_FULL %>%
            filter(`fit$classification`==1)
          TOF_2 <- TOF_FULL %>%
            filter(`fit$classification`==2)
          if(mean(TOF_2$log_TOF) > mean(TOF_1$log_TOF)){
            mylist <- rbind(mylist,TOF_2)
          }else{
            mylist <- rbind(mylist,TOF_1)
          }
        }
      }}}

#Sum plate and then run bioprune and prune outliers from Easysorter package
summedplate <- easysorter::sumplate(mylist,picked_plates = TRUE, directories = FALSE, quantiles = TRUE)

summedplate_1 <- summedplate %>%
  mutate(norm.n=100)
biopruned <- easysorter::bioprune(summedplate_1)

out_pruned <- easysorter::prune_outliers(biopruned)

save(out_pruned, file = "~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Processed_data/20190607pruned_constrain.Rdata")

#Subtract DMSO condition from treatment conditions
biopruned_pr <- biopruned %>%
  dplyr::ungroup() %>%
  dplyr::select(-(date:assay)) %>%
  dplyr::filter(!is.na(condition)) %>%
  tidyr::gather(trait, value, -(plate:col))

control_dr <- biopruned_pr %>%
  dplyr::filter(condition == "DMSO") %>%
  dplyr::group_by(plate,strain, trait) %>%
  dplyr::summarise(control_value = mean(value, na.rm = T))

subtract_dr_to_pc <- biopruned_pr %>%
  dplyr::left_join(., control_dr, by = c("plate","strain", "trait")) %>%
  dplyr::mutate(delta_control = value - control_value) %>%
  dplyr::select(strain, condition, trait, delta_control, plate, row, col)%>%
  dplyr::mutate(u_strain = paste(strain, condition, plate,row,col,sep="_"))%>%
  tidyr::spread(trait, delta_control) %>%
  dplyr::select(-contains("green"),-contains("yellow"), -contains("red"),-contains("cv"),-contains("iqr"), -contains("var"))

save(subtract_dr_to_pc, file = "~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Processed_data/20190607subtract_constrain.Rdata")

