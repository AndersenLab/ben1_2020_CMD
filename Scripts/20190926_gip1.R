library(cegwas)

#get linkage mapped confidence interval from cegwas
x<- snpeff("gip-1")
data.dir <- paste0("~/Dropbox/HTA/Results/2017_anthelmintic_gwa/analysis/ben1_paper/Raw_data/")

#load easysorter data from gwas and filter to Albendazole
load("~/Dropbox/HTA/Results/2017_anthelmintic_gwa/analysis/2018_analysis/Data/2018_all_traits_combined.Rda")

ABZ <- df_merge2 %>%
  select(strain,`Albendazole_q90.EXT_ctrl-regressed`)%>%
  filter(!is.na(`Albendazole_q90.EXT_ctrl-regressed`))

#load in ben-1 variants

ben1_indels <- data.table::fread(paste0(data.dir, "ben1_Indels/20171021_ben1_indels.csv"))%>%
  na.omit()

pr_indels <- ben1_indels %>%
  tidyr::unite(marker, Type, Start, End, sep ="_",remove=F)%>%
  dplyr::filter(comments == "frame shift" | grepl("Exon", location) | Type == "trans")%>%
  dplyr::filter(!grepl("Intron",location,ignore.case = T))%>%
  dplyr::mutate(GT = ifelse(marker == "_NA_NA", "REF", "ALT"))%>%
  dplyr::select(marker, strain = Strain, GT,start = Start,end =End)%>%
  dplyr::distinct(strain, marker, GT,.keep_all=T)

# # # LOAD AND PROCESS SNPS
ben1_snps <- cegwas::snpeff("ben-1",severity = "ALL",elements = "ALL")

ben1_snps_pr <- ben1_snps%>%
  dplyr::select(POS, strain, GT, nt_change, effect, impact)%>%
  tidyr::unite(marker, effect,impact,nt_change,POS,sep="_",remove=F)%>%
  dplyr::select(marker, strain, GT,start = POS)%>%
  dplyr::mutate(end = start)%>%
  dplyr::filter(!grepl("MODIF|syn",marker))

# # # COMBINE SNPS AND INDELS
ben1_variants <- dplyr::bind_rows(ben1_snps_pr,pr_indels)%>%
  dplyr::filter(GT=="ALT")
stop_gain <- x %>%
  filter(GT =="ALT" & impact == "HIGH" & POS == "2361207", FT == "PASS")

#Mark strains with ALT alleles for either gene and both and combinen 
alts_sg_nob1 <- ABZ %>%
  filter(strain %in% stop_gain$strain & !(strain %in% ben1_variants$strain))%>%
  mutate(gt = "stop_gain")
alts_sg_b1 <- ABZ %>%
  filter(strain %in% stop_gain$strain & strain %in% ben1_variants$strain)%>%
  mutate(gt = "B1_SG")
alts_nosg_b1 <- ABZ %>%
  filter(strain %in% ben1_variants$strain & !(strain %in% stop_gain$strain))%>%
  mutate(gt="B1_noSG")
refs_FS <- ABZ %>%
  filter(!(strain %in% stop_gain$strain) & !(strain %in% ben1_variants$strain))%>%
  mutate(gt="ref")

full<- rbind(alts_sg_nob1,alts_sg_b1,alts_nosg_b1,refs_FS)

kruskal.test(full$'Albendazole_q90.EXT_ctrl-regressed', full$gt )

pairwise.wilcox.test(full$'Albendazole_q90.EXT_ctrl-regressed', full$gt)

fit <- aov(full$'Albendazole_q90.EXT_ctrl-regressed' ~ full$gt, data = full, projections = TRUE)
TukeyHSD(fit)

