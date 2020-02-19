#Source the theme 
source("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Scripts/theme.R")
source("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Scripts/scripts.R")
#Pick Colors for strains
cols <- c("N2" = "orange","882"="grey","919"="pink", "920" = "pink", "1081" = "pink", "1075" = "pink", "1325" = "pink","1326"="pink", "1327" = "pink","1328"  = "pink","N2" = "orange", "919" = "pink", "1082" = "pink", "1076" = "pink", "1326" = "pink", "1328" = "pink")
#Parasite strains"
p <- c("N2", "882","919","920", "1325","1326","1327","1328","1081","1082","1076","1075")

##Make plots for R01
#Load ABZ data
load("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Processed_data/20190617subtract_constrain.Rdata")
ABZ_subtract <- subtract_dr_to_pc
rm(subtract_dr_to_pc)

#Load FBZ data
load("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Processed_data/20190618subtract_constrain.Rdata")
FBZ_subtract <- subtract_dr_to_pc
rm(subtract_dr_to_pc)

##Figure 1 q90.EXT

load("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Processed_data/20190617subtract_pruned.Rdata")

r01_strains <- c("N2","882")

abz_r01 <- x%>%
  dplyr::filter(strain %in% r01_strains, trait == "q90.EXT")%>%
  dplyr::filter(!is.na(trait))
rm(x)
abz_norm <- abz_r01%>%
  dplyr::select(strain,condition, trait,phenotype)%>%
  dplyr::mutate(q90EXT_norm = ((abz_r01$phenotype - min(abz_r01$phenotype)) / (max(abz_r01$phenotype) - min(abz_r01$phenotype))))

abz_plot <- abz_norm%>%
  filter(condition!="200",
         !is.na(strain))%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "882","919","920","1081","1082","1327","1328", "1325","1326","1075","1076","1137","1138","1139","1140","1314","1317","1097","1098")))%>%
  dplyr::mutate(con_fix = factor(condition, levels = c("DMSO","6_25","12_5","25","50","100","200")))%>%
  ggplot()+
  aes(x=con_fix, y = q90EXT_norm)+
  geom_smooth(alpha=0.1, se=TRUE, aes(group=fancy_strain, colour=fancy_strain), method = 'auto')+
  theme_cowplot(24)+
  ylim(0.25,0.95)+
  ylab("ABZ response")+
  scale_x_discrete(labels=c("DMSO" = "DMSO", "6_25" = "6.25 µM","12_5" = "12.5 µM","25"="25 µM","50"="50 µM","100"="100 µM"))+
  scale_fill_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y","1138"="A185P","1140"="A185P","1138"="Q131L","1098"="S145F", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V","1137"="Q131L","1139"="A185P","1314"="M257I","1317"="M257I","1097"="S145F"), values = cols)+
  scale_color_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y","1138"="A185P","1140"="A185P","1138"="Q131L","1098"="S145F", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V","1137"="Q131L","1139"="A185P","1314"="M257I","1317"="M257I","1097"="S145F"), values = cols)+
  base_theme+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(face = "bold"),
        axis.text.y=element_text(face = "bold"),
        axis.title.y=element_text(face = "bold"),
        legend.key.size = unit(0.4,"in"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
        legend.position = "None")

load("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Processed_data/20190618subtract_pruned.Rdata")

r01_strains <- c("N2","882")

fbz_r01 <- x%>%
  dplyr::filter(strain %in% r01_strains, trait == "q90.EXT")%>%
  dplyr::filter(!is.na(trait))
rm(x)
fbz_norm <- fbz_r01%>%
  dplyr::select(strain,condition, trait,phenotype)%>%
  dplyr::mutate(q90EXT_norm = ((fbz_r01$phenotype - min(fbz_r01$phenotype)) / (max(fbz_r01$phenotype) - min(fbz_r01$phenotype))))

fbz_plot <- fbz_norm%>%
  filter(condition!="200",
         !is.na(strain))%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "882","919","920","1081","1082","1327","1328", "1325","1326","1075","1076","1137","1138","1139","1140","1314","1317","1097","1098")))%>%
  dplyr::mutate(con_fix = factor(condition, levels = c("DMSO","6_25","12_5","25","50","100","200")))%>%
  ggplot()+
  aes(x=con_fix, y = q90EXT_norm)+
  geom_smooth(alpha=0.1, se=TRUE, aes(group=fancy_strain, colour=fancy_strain), method = 'auto')+
  theme_cowplot(24)+
  ylim(0.25,0.95)+
  ylab("FBZ response")+
  scale_x_discrete(labels=c("DMSO" = "DMSO", "6_25" = "6.25 µM","12_5" = "12.5 µM","25"="25 µM","50"="50 µM","100"="100 µM"))+
  scale_fill_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y","1138"="A185P","1140"="A185P","1138"="Q131L","1098"="S145F", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V","1137"="Q131L","1139"="A185P","1314"="M257I","1317"="M257I","1097"="S145F"), values = cols)+
  scale_color_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y","1138"="A185P","1140"="A185P","1138"="Q131L","1098"="S145F", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V","1137"="Q131L","1139"="A185P","1314"="M257I","1317"="M257I","1097"="S145F"), values = cols)+
  base_theme+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(face = "bold"),
        axis.text.y=element_text(face = "bold"),
        axis.title.y=element_text(face = "bold"),
        legend.key.size = unit(0.4,"in"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
        legend.position = "None")

cowplot::plot_grid(abz_plot, fbz_plot, labels = c('A', 'B'), label_size = 12, nrow = 2,ncol=1)
ggsave("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/plots/R01_figures/Figure1.png", plot = last_plot(),width = 5,height = 4, units = "in")
  


##Figure 2 Mean TOF

load("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Processed_data/20190603para_regressed.Rdata" )
abz_V2 <- para_regressed
rm(para_regressed)

abz_v2_plot <- abz_V2 %>%
  filter(trait=="mean.TOF",
         !is.na(strain))%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "882","1075","1076","1081","1082","1325","1326","1327","1328","919","920")))%>%
  ggplot()+
  aes(x=fancy_strain, y = phenotype)+
  geom_jitter(width = 0.1)+
  geom_boxplot(aes(fill=fancy_strain, alpha=0.1),outlier.shape = NA)+
  theme_cowplot(24)+
  ylab("Animal Length (ABZ)")+
  scale_x_discrete(labels=c("N2" = "WT", "882" = "Del","919" = "F200Y","920"="F200Y","1081"="E198A","1076"="F167Y", "1325" = "E198L", "1327" = "E198V", "1082" = "E198A","1075"="F167Y"))+
  scale_fill_manual(name = "Strain", labels = c("N2" = "N2","882"="Del","919"="F200Y","1140"="A185P","1138"="Q131L","1098"="S145F", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V","1137"="Q131L","1314"="Q131L","1317"="M257I","1097"="S145F"), values = cols)+
  base_theme+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(face = "bold"),
        axis.text.y=element_text(face = "bold"),
        axis.title.y=element_text(face = "bold"),
        legend.key.size = unit(0.4,"in"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
        legend.position = "None")

load("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Processed_data/20190606para_regressed.Rdata")
fbz_V3 <- para_regressed

fbz_v3_plot <- fbz_V3 %>%
  filter(trait=="mean.TOF",
         !is.na(strain))%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "882","1075","1076","1081","1082","1325","1326","1327","1328","919","920")))%>%
  ggplot()+
  aes(x=fancy_strain, y = phenotype)+
  geom_jitter(width = 0.1)+
  geom_boxplot(aes(fill=fancy_strain,alpha=0.1),outlier.shape = NA)+
  theme_cowplot(24)+
  ylab("Animal Length (FBZ)")+
  scale_x_discrete(labels=c("N2" = "WT", "882" = "Del","919" = "F200Y","920"="F200Y","1139"="A185P","1081"="E198A","1138"="Q131L","1076"="F167Y", "1325" = "E198L", "1327" = "E198V", "1082" = "E198A","1075"="F167Y","1137"="Q131L","1097"="S145F","1098"="S145F","1140" = "A185P","1317"="M257I","1328" = "E198V"))+
  scale_fill_manual(name = "Strain", labels = c("N2" = "N2","882"="Del","919"="F200Y","1140"="A185P","1138"="Q131L","1098"="S145F", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V","1137"="Q131L","1314"="Q131L","1317"="M257I","1097"="S145F"), values = cols)+
  base_theme+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(face = "bold"),
        axis.text.y=element_text(face = "bold"),
        axis.title.y=element_text(face = "bold"),
        legend.key.size = unit(0.4,"in"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
        legend.position = "None")
cowplot::plot_grid(abz_v2_plot, fbz_v3_plot, labels = c('A', 'B'), label_size = 12, nrow = 2,ncol=1)
ggsave("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/plots/R01_figures/Figure2.png", plot = last_plot(),width = 7.5,height = 5, units = "in")

#Figure 3 median.TOF
load("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Processed_data/Rescue_regressed.Rdata")
colsgfp <- c("ECA884" = "grey", "1108_noGFP" = "grey", "1108_GFP" = "green", "N2" = "orange")

rescue_data %>%
  filter(trait=="median.TOF",
         !is.na(strain)) %>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "ECA884", "1108_GFP", "1108_noGFP")))%>%
  ggplot()+
  aes(x=fancy_strain, y = phenotype)+
  geom_jitter(size=1,width = 0.25)+
  geom_boxplot(aes(fill=fancy_strain, alpha=0.1),outlier.shape = NA)+
  ylab("Residual Animal Length (ABZ)")+
  theme_cowplot(24)+
  scale_x_discrete(labels=c("ECA884" = "Del", "1108_noGFP" = "Rescue(-)","1108_GFP" = "Rescue(+)", "N2" = "WT"))+
  scale_fill_manual(name = "Strain", labels = c("ECA884" = "∆ben-1", "ECA1108_noGFP" = "Rescue GFP(-)", "ECA1108_GFP" = "Rescue GFP(+)", "N2"), values = colsgfp)+
  base_theme+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(face = "bold",size=10,angle = 45,vjust=0.7),
        axis.text.y=element_text(face = "bold",hjust=1, ),
        axis.title.y=element_text(face = "bold", vjust=-1.5),
        legend.key.size = unit(0.4,"in"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position = "None")
ggsave("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/plots/R01_figures/Figure3.png",plot = last_plot(),width = 3,height = 3, units = "in")


##Figure 5 Thiabendazole ben-1 regressed (mean.TOF), Albendazole ben-1 removed q90EXT, Albendazole ben-1 regressed q90TOF

source("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Scripts/scripts.R")

df_TBZmeanTOF <- rio::import("~/Dropbox/HTA/Results/2017_anthelmintic_gwa/analysis/2018_analysis/Data/Mapping_results/ben-1_regressed/mapping_Rdata/Thiabendazole-mean-TOF-mapping.Rda")

TBZmeanTOF<-plot_manhatten(df_TBZmeanTOF, option = "mid")

df_ABZq90 <- rio::import("~/Dropbox/AndersenLab/LabFolders/Clay/ben1strains_removed_mapping/Mappings/Data/Albendazole_q90.EXT_processed_mapping.tsv")

ABZq90EXT <- plot_manhatten(df_ABZq90, option = "bot")

df_ABZq90TOF <- rio::import("~/Dropbox/HTA/Results/2017_anthelmintic_gwa/analysis/2018_analysis/Data/Mapping_results/ben-1_regressed/mapping_Rdata/Albendazole-q90-TOF-mapping.Rda")

ABZq90TOF <- plot_manhatten(df_ABZq90TOF, option = "top")

manhatten_plots <- plot_grid(ABZq90TOF,TBZmeanTOF,ABZq90EXT,ncol = 1, nrow = 3, labels = c("A","B","C"))
ggsave("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/plots/R01_figures/Figure5.png", plot = manhatten_plots, width = 7.5, height = 6.5, units = "in")


#Figure 6 mean.TOF
load("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Processed_data/20190923_125.Rdata")
load("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Processed_data/20190924_125.Rdata")
load("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Processed_data/20190923_625.Rdata")
load("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Processed_data/20190924_625.Rdata")

bleach_625_d12 <- rbind(bleach_625_d1,bleach_625_d2)

assay_regressed <- regress(bleach_625_d12, assay = TRUE)

full_thia <- regress(assay_regressed)

traita <- "q25.TOF"
cols <- c("N2" = "orange", "CB4856" = "blue", "238" = "grey", "239"="grey")
thia <- full_thia%>%
  filter(trait == traita,
         !is.na(condition),
         !(grepl("B", condition)),
         !is.na(strain))%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "CB4856", "238", "239")))%>%
  dplyr::mutate(con_fix = factor(condition, levels = c("DMSO","31_25","62_5","125","250")))%>%
  ggplot()+
  aes(x=fancy_strain, y = phenotype)+
  geom_jitter(size=1,width = 0.1)+
  geom_boxplot(aes(fill=strain),alpha=0.2,outlier.shape = NA)+
  theme_cowplot(24)+
  ylab("Animal Length [TBZ]")+
  scale_x_discrete(labels=c("N2" = "N2", "CB4856" = "CB4856","238" = "N2[CB4856 V]","239"="CB4856[N2 V]"))+
  scale_fill_manual(name = "Strain", labels = c("N2" = "WT", "CB4856" = "CB4856", "238" = "238", "239" = "239"), values = cols)+
  base_theme+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(face = "bold",size=8,angle = 45,vjust=0.7),
        axis.text.y=element_text(face = "bold",hjust=1, ),
        axis.title.y=element_text(face = "bold", vjust=-1.5),
        legend.key.size = unit(0.4,"in"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position = "None")
ggsave("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/plots/R01_figures/Figure6.png", plot = thia, device = "png",width = 3,height = 3, units = "in")

##Figure 7

kaplan <- read_excel("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Processed_data/kaplan_EC50.xlsx")

kaplan_adj <-kaplan %>%
  mutate(farm = "Farm")%>%
  mutate(adjusted = as.numeric(kaplan$'PRISM-BZ_Conc'))

kplan_plot <- ggplot(kaplan_adj)+
  aes(x= farm, y=adjusted)+
  geom_quasirandom(cex = 0.1)+
  scale_y_log10()+
  geom_hline(yintercept = 0.117, linetype = 1, color="grey")+
  geom_hline(yintercept = 0.234, linetype = 1, color ="grey")+
  labs(y = expression(paste(bold("EC"[50]*" (µM)"))))+
  base_theme+
  theme_cowplot(24)+
  theme(
    axis.title.x=element_blank(),
    axis.line = element_line(size = 1),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y=element_text(face = "bold",hjust=1,size=8 ),
    axis.title.y=element_text(face = "bold", vjust=-1.5,size = 8),
    legend.key.size = unit(0.4,"in"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    legend.position = "None")

ggsave("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/plots/R01_figures/Figure7.png", plot = kplan_plot, device = "png", width =1.5, height = 2, units="in")

#Figure 8  
load("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Processed_data/Lifestagebetas.Rdata")
ancestry.colours <- c('#F2F3F4', '#222222', '#F3C300', '#875692', '#F38400', '#A1CAF1', '#BE0032', '#C2B280', '#848482', '#008856', '#E68FAC', '#0067A5', '#F99379', '#604E97', '#F6A600', '#B3446C', '#DCD300', '#882D17', '#8DB600', '#654522', '#E25822', '#2B3D26')
lifestage <- ggplot(data=BetaTub_LS, aes(x=Life.Stage, y=FPKM.value, col = id)) +
  geom_line(size=1) +
  geom_errorbar(aes(ymin=FPKM.value-sd, ymax=FPKM.value+sd), width=.1,size=0.5) +
  scale_colour_manual(values = ancestry.colours[2:10], name = "Beta Tublin") +
  ylab("Expression Levels [FPKM]") +
  base_theme+
  theme_cowplot(24)+
  xlab("Life Stage") +
  scale_x_continuous(breaks=c(1:6), labels=c("Embryo" , "L1", "L2" , "L3" , "L4", "Adult"))+
  theme(
    axis.title.x=element_blank(),
    axis.line = element_line(size = 1),
    axis.text.x = element_text(face = "bold",size =5),
    axis.text.y=element_text(face = "bold",hjust=1 , size=6),
    axis.title.y=element_text(face = "bold",vjust=-0.5, size=6),
    legend.key.size = unit(0.4,"in"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    legend.position = "none")
library(monocle)
library(tidyverse)
library(ggplot2)

load("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Processed_data/Cao_et_al_2017_vignette.RData")
x <-show.expr.info("ben-1", "tissue")

ben1_tissue <- x %>%
  ggplot()+
  aes(x=facet,y=tpm/100)+
  geom_col()+
  base_theme+
  theme_cowplot(24)+
  labs(y = "Expression Levels [TPM/1e3]")+
  scale_x_discrete(labels = c("Neurons" = "Neurons", "Gonad" = "Gonad","Glia"="Glia","Hypodermis"="Hypodermis","Intestine"="Intestine","Pharynx"="Pharynx","Body wall muscle" = "Muscle"))+
  theme(axis.title.x=element_blank(),
        axis.line = element_line(size = 1),
        axis.text.x = element_text(face = "bold", size = 5),
        axis.text.y=element_text(face = "bold",hjust=1,size = 6),
        axis.title.y=element_text(face = "bold",vjust=-0.5, size=6),
        legend.key.size = unit(0.4,"in"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))
expression <- plot_grid(lifestage,ben1_tissue,nrow = 1)
ggsave("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/plots/R01_figures/Figure8.png",plot = expression, width = 6,height = 1.5,units = "in")
