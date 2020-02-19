library(cegwas2)

data.dir <- paste0("~/Dropbox/HTA/Results/2017_anthelmintic_gwa/analysis/ben1_paper/Raw_data/")

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
#get linkage mapped confidence interval from cegwas
x<- snpeff("V:15734606-16165935")

#Get strains with CB allele

ALT_strains_cyp_CB <- x %>%
  filter(GT=="ALT" & gene_name=="cyp-35D1" & aa_change=="p.Lys267Glu")%>%
  distinct(strain)

#load easysorter data from gwas and filter to thiabendazole
load("~/Dropbox/HTA/Results/2017_anthelmintic_gwa/analysis/2018_analysis/Data/2018_control_regressed_data_outliers_removed.Rda")

thia <- ctrl_regressed %>%
  select(strain,Thiabendazole_q90.TOF)%>%
  filter(!is.na(Thiabendazole_q90.TOF))

#mark strains with CB allele and ben-1 allele
alts_cyp_CB <- thia %>%
  filter(strain %in% ALT_strains_cyp_CB$strain )%>%
  mutate(gt = "CB4856 Allele")

refs_cyp <- thia%>%
  filter(!(strain %in% ALT_strains_cyp_CB$strain))%>%
  mutate(gt ="Other Alleles")


al_re_cyp <- rbind(refs_cyp,alts_cyp_CB)

colors_grey <- c("CB4856 Allele" = "blue", "CB4856 Allele with additional changes" = "grey","CBVar" = "blue", "Other Alleles" = "grey","B1var_CBvar"="blue", "REF" = "grey","B1var" = "grey","novar"="grey")
colors_cyp <- c("CB4856 Allele" = "blue", "CB4856 Allele with additional changes" = "darkgreen","CBVar" = "blue", "Other Alleles" = "grey28","B1var_CBvar"="green", "REF" = "grey","B1var" = "red","novar"="grey")

b1_var <- al_re_cyp %>%
  filter(strain %in% ben1_variants$strain) %>%
  filter(gt == "CB4856 Allele")

ben_1_var_lab <- al_re_cyp %>%
  mutate(ben_1 = case_when(strain %in% b1_var$strain ~ "B1var_CBvar",
                           strain %in% ben1_variants$strain ~ "B1var",
                           gt == "CB4856 Allele" ~ "CBVar",
                           !(strain %in% ben1_variants$strain) ~ "novar"))

kruskal.test(ben_1_var_lab$Thiabendazole_q90.TOF, ben_1_var_lab$ben_1)

pairwise.wilcox.test(ben_1_var_lab$Thiabendazole_q90.TOF, ben_1_var_lab$ben_1)


