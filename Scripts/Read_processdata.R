library(easysorter)
library(COPASutils)
library(dplyr)
library(tidyverse)
library(cowplot)
library(mclust)
parafull <- c("N2","882","919","920","1075","1076","1081","1082","1325","1326","1327","1328")
traita <- c("median.EXT")
cols <- c("N2" = "orange","882"="grey","919"="red","920"="red","1325"="yellow","1326"="yellow","1327"="purple","1328"="purple","1082" = "blue","1081" = "blue", "1076" = "green","1075" = "green")

###Read in ABZ dose data from 20190617 
dir <- c("~/Desktop/ben1_2020_CMD/Raw_data/S1/")

raw <- easysorter::read_data(dir)

raw_noncontam <- easysorter::remove_contamination(raw)

summedplate <- easysorter::sumplate(raw_noncontam, directories = FALSE, quantiles = TRUE,v3_assay = TRUE)

biopruned <- easysorter::bioprune(summedplate)

out_pruned <- easysorter::prune_outliers(biopruned)%>%
  dplyr::filter(strain %in% parafull)%>%
  dplyr::filter(trait == traita)

DMSO_pheno <- out_pruned %>%
  dplyr::filter(condition == "DMSO") %>%
  dplyr::group_by(strain,date) %>%
  dplyr::summarise(cntrl_pheno = mean(phenotype, na.rm = T))%>%
  dplyr::ungroup()%>%
  dplyr::select(strain,cntrl_pheno)

subtracted_dose <- out_pruned %>%
  dplyr::left_join(., DMSO_pheno, by = c("strain"))%>%
  dplyr::mutate(sub_cntrl = phenotype - cntrl_pheno)%>%
  dplyr::mutate(perc_DMSO = phenotype/cntrl_pheno)

save(subtracted_dose,file = "~/Desktop/ben1_2020_CMD/Processed_data/ABZv3_dose.RData")

rm(dir,biopruned,DMSO_pheno,out_pruned,raw,raw_noncontam,subtracted_dose,summedplate)

###Read and process ABZ high-replication canonical alleles

dirs_abzH <- list.files("S3",full.names = TRUE)
raw_abzH <- easysorter::read_data(dirs_abzH)

raw_noncontam_abzH <- easysorter::remove_contamination(raw_abzH)


summedplate_abzH <- easysorter::sumplate(raw_noncontam_abzH, directories = TRUE, quantiles = TRUE,v3_assay = TRUE)

biopruned_abzH <- easysorter::bioprune(summedplate_abzH) 
save(biopruned_abzH, file = "~/Desktop/ben1_2020_CMD/Processed_data/DMSOABZ.RData")
assay_regressed_abzH <- easysorter::regress(biopruned_abzH, assay = TRUE)

out_pruned_abzH <- easysorter::prune_outliers(assay_regressed_abzH)

regressedabz_abzH <- out_pruned_abzH %>%
  dplyr::filter(strain %in% c("N2","882","919", "1076", "1082"))%>%
  dplyr::filter(trait == traita)%>%
  dplyr::mutate(assay = "A")%>%
  easysorter::regress()

save(regressedabz_abzH,file = "~/Desktop/ben1_2020_CMD/Processed_data/ABZv3_highrepO.RData")

###Read in FBZ dose from 20190618

dir_fbz <- c("~/Desktop/ben1_2020_CMD/Raw_data/S2/")

raw_fbz <- easysorter::read_data(dir_fbz)

raw_noncontam_fbz <- easysorter::remove_contamination(raw_fbz)

summedplate_fbz <- easysorter::sumplate(raw_noncontam_fbz, directories = FALSE, quantiles = TRUE,v3_assay = TRUE)

biopruned_fbz <- easysorter::bioprune(summedplate_fbz)

out_pruned_fbz <- easysorter::prune_outliers(biopruned_fbz)%>%
  dplyr::filter(strain %in% parafull)%>%
  dplyr::filter(trait == traita)

DMSO_pheno_fbz <- out_pruned_fbz %>%
  dplyr::filter(condition == "DMSO") %>%
  dplyr::group_by(strain,date) %>%
  dplyr::summarise(cntrl_pheno = mean(phenotype, na.rm = T))%>%
  dplyr::ungroup()%>%
  dplyr::select(strain,cntrl_pheno)

subtracted_dose_fbz <- out_pruned_fbz %>%
  dplyr::left_join(., DMSO_pheno_fbz, by = c("strain"))%>%
  dplyr::mutate(sub_cntrl = phenotype - cntrl_pheno)%>%
  dplyr::mutate(perc_DMSO = phenotype/cntrl_pheno)

save(subtracted_dose_fbz,file = "~/Desktop/ben1_2020_CMD/Processed_data/FBZv3_dose.RData")

rm(dir_fbz,biopruned_fbz,DMSO_pheno_fbz,out_pruned_fbz,raw_fbz,raw_noncontam_fbz,subtracted_dose_fbz,summedplate_fbz)

### Read and process FBZ high-replication old alleles

dir_fbzH <- c("~/Desktop/ben1_2020_CMD/Raw_data/S4/")

raw_fbzH <- easysorter::read_data(dir_fbzH)

raw_noncontam_fbzH <- easysorter::remove_contamination(raw_fbzH)

summedplate_fbzH <- easysorter::sumplate(raw_noncontam_fbzH, directories = FALSE, quantiles = TRUE, v3_assay = TRUE)

biopruned_fbzH <- easysorter::bioprune(summedplate_fbzH)

regressed_fbzH <- biopruned_fbzH %>%
  easysorter::prune_outliers(., iqr = FALSE)%>%
  dplyr::filter(strain %in% c("N2","882","920","1076","1082"))%>%
  easysorter::regress(.,assay = FALSE)%>%
  dplyr::filter(trait == traita)

save(regressed_fbzH, file = "~/Desktop/ben1_2020_CMD/Processed_data/FBZv3_highrepO.RData")

### Read and process ABZ high-replication new alleles 
 
dirs_abzHN <- list.files("S3", full.names = TRUE)
raw_abzHN <- easysorter::read_data(dirs_abzHN)


raw_noncontam_abzHN <- easysorter::remove_contamination(raw_abzHN)

summedplate_abzHN <- easysorter::sumplate(raw_noncontam_abzHN, directories = TRUE, quantiles = TRUE,v3_assay = TRUE)

biopruned_abzHN <- easysorter::bioprune(summedplate_abzHN) 

assay_regressed_abzHN <- easysorter::regress(biopruned_abzHN, assay = TRUE)

out_pruned_abzHN <- easysorter::prune_outliers(assay_regressed_abzHN)

regressedabz_abzHN <- out_pruned_abzHN %>%
  dplyr::filter(strain %in% c("N2","882","1325", "1327"))%>%
  dplyr::filter(trait == traita)%>%
  dplyr::mutate(assay = "A")%>%
  easysorter::regress()

save(regressedabz_abzHN,file = "~/Desktop/ben1_2020_CMD/Processed_data/ABZv3_highrepN.RData")

### Read and process FBZ high-replication new alleles 

dir_fbzHN <- c("~/Desktop/ben1_2020_CMD/Raw_data/S4/")

raw_fbzHN <- easysorter::read_data(dir_fbzHN)

raw_noncontam_fbzHN <- easysorter::remove_contamination(raw_fbzHN)

summedplate_fbzHN <- easysorter::sumplate(raw_noncontam_fbzHN, directories = FALSE, quantiles = TRUE, v3_assay = TRUE)

biopruned_fbzHN <- easysorter::bioprune(summedplate_fbzHN)

regressed_fbzHN <- biopruned_fbzHN %>%
  easysorter::prune_outliers(., iqr = FALSE)%>%
  dplyr::filter(strain %in% c("N2","882","1325","1327"))%>%
  easysorter::regress(.,assay = FALSE)%>%
  dplyr::filter(trait == traita)

save(regressed_fbzHN, file = "~/Desktop/ben1_2020_CMD/Processed_data/FBZv3_highrepN.RData")

### Calculate EC50

load("~/Desktop/ben1_2020_CMD/Processed_data/FBZv3_dose.RData")
plate_1strains <- c("N2","882","919","920","1081","1082","1075","1076","1325","1326","1327","1328")
fixcon <- subtracted_dose_fbz%>%
  dplyr::filter(!(is.na(strain)))%>%
  dplyr::filter(condition != "200")%>%
  dplyr::filter(!(strain == "1328" & plate == 14))%>%  #Filtered because no DMSO condition for this plate
  dplyr::filter(!(strain == "1328" & plate == 19))%>% #Filtered because concentrations above 12.5µM not included
  dplyr::mutate(cond = case_when(condition == "DMSO" ~ 0,
                                 condition == "6_25" ~ 6.25,
                                 condition == "12_5" ~ 12.5,
                                 condition == "25" ~ 25,
                                 condition == "50" ~ 50,
                                 condition == "100" ~ 100))%>%
  dplyr::filter(trait == "median.EXT")%>%
  dplyr::filter(strain %in% plate_1strains)%>%
  dplyr::mutate(perc_DMSO = phenotype / cntrl_pheno)

collectEC <- data.frame(strain = character(),
                        EC50 = numeric(),
                        plate = numeric(),
                        par_b = numeric(),
                        par_c = numeric(),
                        model = character(),
                        stringsAsFactors = FALSE)

for (i in plate_1strains){
  strain_1 <- fixcon %>%
    dplyr::filter(strain == i)
  for (j in (1:20)) {
    plates <- strain_1 %>%
      dplyr::filter(plate == j)
    skip_to_next <- FALSE
    tryCatch(all <- drc::drm(phenotype ~ cond, strain, 
                             data = plates,
                             fct = LL.3()), error = function(e) { skip_to_next <<- TRUE})
    if(skip_to_next) { 
      collectEC <- add_row(collectEC, 
                           strain = i,
                           EC50 = NA,
                           par_b = NA,
                           par_c = NA,
                           model= "LL.3",
                           plate = j)
      next } 
    tryCatch( ec50N <- ED(all,50), error = function(e) { skip_to_next <<- TRUE})
    Coef <- all$coefficients
    collectEC <- add_row(collectEC, 
                         strain = i,
                         EC50 = ec50N[1,1],
                         par_b = Coef[1],
                         par_c = Coef[2],
                         model= "LL.3",
                         plate = j)
  }
}
save(collectEC, file = "~/Desktop/ben1_2020_CMD/Processed_data/EC50.RData")

#Save regressed high-throughput FBZ for supplemental figure 3
regressed_fbzHA <- biopruned_fbzH%>%
  dplyr::filter(strain %in% plate_1strains)%>%
  dplyr::ungroup()%>%
  dplyr::mutate(assay = case_when(plate %in% c(1,2,3,4,5,6,7,8,9,10) ~ "A",
                                  plate %in% c(11,12,13,14,15,16,17,18,19,20) ~ "B"))%>%
  easysorter::regress(., assay = TRUE)%>%
  easysorter::prune_outliers(.,iqr=FALSE)%>%
  easysorter::regress(.,assay = FALSE)%>%
  dplyr::filter(trait == traita)
save(regressed_fbzHA, file = "~/Desktop/ben1_2020_CMD/Processed_data/FBZv3_highrepA.RData")

