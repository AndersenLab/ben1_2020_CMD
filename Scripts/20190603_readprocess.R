library(easysorter)
library(COPASutils)
library(dplyr)
library(tidyverse)
library(cowplot)

###Read in day 1 and analyze
setwd("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Processed_data/")
dir <- c("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Raw_data/20190603_ABZ/")

raw <- read_data(dir)

raw_noncontam <- remove_contamination(raw)

summedplate <- sumplate(raw_noncontam, directories= FALSE, quantiles = TRUE)

biopruned <- bioprune(summedplate)

out_pruned <- prune_outliers(biopruned)

full <- regress(out_pruned)

p <- c("N2", "882","919","1082","1075")
plate_1 <- c("1","2","3","4","5","6","7","8","9","10")
para_full <- out_pruned %>%
  dplyr::filter(strain %in% p)%>%
  dplyr::filter(plate %in% plate_1)

para_regressed <- regress(para_full)

save(para_regressed, file = "~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Processed_data/20190603para_regressed.Rdata")

save(out_pruned, file = "~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Processed_data/20190603pruned.Rdata")

save(full, file = "~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Processed_data/20190603regressed.Rdata")


