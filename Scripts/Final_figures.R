#Source the theme and required libraries
library(tidyverse)
library(cowplot)
library(sjPlot)
library(drc)
library(ggpubr)
library(glue)
source("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Scripts/theme.R")
parafull <- c("N2","882","919","920","1075","1076","1081","1082","1325","1326","1327","1328")
para_con <- c("N2","882","919","1076","1082")
para_new <- c("N2","882", "1325","1327")
setwd("~/Desktop/ben1_2020_CMD/data/")

#Choose colors for strains 
cols <- c("N2" = "orange","882"="grey","919"="red","920"="red","1325"="yellow","1326"="yellow","1327"="purple","1328"="purple","1082" = "blue","1081" = "blue", "1076" = "green","1075" = "green")
para <- c("N2","882","919","1082","1076","1325","1327")
plate_1strains <- c("N2","882","919","920","1081","1082","1075","1076","1325","1326","1327","1328")
traita <- c("median.EXT")

##Figure 1
#Figure 1 A

load("S1.RData")
abz_outpruned <- subtracted_dose

abz_outprunedmut <- abz_outpruned%>%
  dplyr::filter(condition != "200")%>%
  dplyr::filter(trait == traita)%>%
  dplyr::filter(strain %in% para_con)%>%
  dplyr::mutate(numcon = case_when(condition == "DMSO" ~ 0,
                                   condition == "6_25" ~ 6.25,
                                   condition == "12_5" ~ 12.5,
                                   condition == "25" ~ 25,
                                   condition == "50" ~ 50,
                                   condition == "100" ~ 100,
                                   condition == "200" ~ 200))%>%
  dplyr::group_by(strain, condition)%>%
  dplyr::mutate(meancondition = mean(sub_cntrl))%>%
  dplyr::mutate(SND = sd(sub_cntrl))

abzmedian.EXT_subtractedold <- abz_outprunedmut%>%
  dplyr::filter(!is.na(strain))%>%
  dplyr::filter(strain %in% c("N2","919","920","882","1081","1082", "1075","1076","1325","1326","1327","1328"))%>%
  dplyr::filter(condition != "200")%>%
  dplyr::filter(numcon %in% c(0,6.25,12.5,25,50,100,200))%>%
  dplyr::filter(trait == traita)%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "882","919","920","1081","1082","1327","1328", "1325","1326","1075","1076")))%>%
  dplyr::mutate(con_fix = factor(condition, levels = c("DMSO","6_25","12_5","25","50","100","200")))%>%
  ggplot()+
  aes(x=numcon, y = meancondition)+
  geom_errorbar(aes(color = fancy_strain, ymin=meancondition - SND, ymax=meancondition + SND), width=2, size = 1)+
  geom_line(aes(x = numcon, y = meancondition, colour= fancy_strain), size = 1)+
  theme_cowplot(12)+
  scale_y_continuous(limits = c(-600,150))+
  ylab("Difference in OD")+
  xlab("Concentration (µM)")+
  scale_x_continuous(breaks = c(0,6.25,12.5,25,50,100),labels = c("0","6.25","12.5","25","50","100"),expand = c(0,0.5))+
  scale_fill_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V"), values = cols)+
  scale_color_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V"), values = cols)+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(size=12,face = "bold",angle = 90,vjust = 0.5,hjust=0.5,margin = unit(c(0,0,0,0),units = "in")),
        axis.text.y = element_text(size=12, face = "bold",margin = unit(c(0 ,0,0,0),units = "in")),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size = 0.5),
        axis.title.y=element_text(size=12,face = "bold"),
        legend.key.size = unit(1,"in"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        plot.margin = unit(c(0, 0, 0, 0), "in"),
        legend.position = "None")
#Figure 1 B

load("S2.RData")
abz_stats <- aov(phenotype ~ strain, data = regressedabz_abzH)%>%
  rstatix::tukey_hsd()%>%
  dplyr::mutate(p.adj.signif = case_when(p.adj > 0.05 ~ "ns",
                                         (p.adj < 0.05 & p.adj > 0.01) ~"*",
                                         (p.adj < 0.01 & p.adj > 0.001) ~"**",
                                         (p.adj < 0.001 & p.adj > 0.0001) ~"***",
                                         (p.adj < 0.0001 & p.adj > 0.00001) ~"****",
                                         (p.adj < 0.00001) ~ "*****"))%>%
  dplyr::filter(grepl("N2",group1) | grepl("N2",group2))


regressedabz_medianEXTold <-regressedabz_abzH%>%
  dplyr::filter(!is.na(strain))%>%
  dplyr::filter(trait == traita)%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "882","919","920","1081","1082","1075","1076","1327","1328", "1325","1326")))%>%
  ggplot()+
  aes(x=fancy_strain, y = phenotype)+
  geom_jitter(width = 0.1,size=0.1)+
  geom_boxplot(aes(fill=strain, alpha = 0.1), outlier.shape = NA)+
  theme_cowplot(8)+
  scale_y_continuous(limits = c(-375, 225), sec.axis = dup_axis(name = "Albendazole"))+
  stat_pvalue_manual(abz_stats, label = "p.adj.signif",xmax="group1", y.position = c(150),remove.bracket = TRUE,size = 4)+
  ylab(glue("Normalized OD"))+
  scale_x_discrete(labels=c("N2" = "N2", "882" = "Deletion","919" = "F200Y","1082"="E198A","1327"="E198V","1325"="E198L","1076"="F167Y"))+
  scale_fill_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V"), values = cols)+
  scale_color_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V"), values = cols)+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=12, face = "bold",margin = unit(c(0 ,0,0,0),units = "in")),
        axis.title.y=element_text(size = 12,face = "bold",margin = unit(c(0,0.075,0,0),units = "in")),
        legend.key.size = unit(0.4,"in"),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(margin = unit(c(0,0,0,0),units = "in")),
        axis.ticks.y.right = element_blank(),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        plot.margin = unit(c(0.1, 0.1, 0, 0), "in"),
        legend.position = "None")

#Figure 1 C

load("S4.RData")
fbz_outpruned <- subtracted_dose_fbz

fbz_outprunedmut <- fbz_outpruned%>%
  dplyr::filter(condition != "200")%>%
  dplyr::filter(trait == traita)%>%
  dplyr::mutate(numcon = case_when(condition == "DMSO" ~ 0,
                                   condition == "6_25" ~ 6.25,
                                   condition == "12_5" ~ 12.5,
                                   condition == "25" ~ 25,
                                   condition == "50" ~ 50,
                                   condition == "100" ~ 100))%>%
  dplyr::group_by(strain, condition)%>%
  dplyr::mutate(meancondition = mean(sub_cntrl))%>%
  dplyr::mutate(SND = sd(sub_cntrl))

FBZmedian.EXT_subtractedold <-fbz_outprunedmut%>%
  dplyr::filter(!is.na(strain))%>%
  dplyr::filter(strain %in% c("N2","919","882","1082", "1076"))%>%
  dplyr::filter(condition != "200")%>%
  dplyr::filter(numcon %in% c(0,6.25,12.5,25,50,100))%>%
  dplyr::filter(trait == traita)%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "882","919","920","1081","1082","1327","1328", "1325","1326","1075","1076")))%>%
  dplyr::mutate(con_fix = factor(condition, levels = c("DMSO","6_25","12_5","25","50","100")))%>%
  ggplot()+
  aes(x=numcon, y = meancondition)+
  geom_errorbar(aes(color = fancy_strain, ymin=meancondition - SND, ymax=meancondition + SND), width=2, size = 1)+
  geom_line(aes(x = numcon, y = meancondition, colour= fancy_strain), size = 1)+
  theme_cowplot(12)+
  scale_y_continuous(limits = c(-750,100))+
  ylab("Difference in OD")+
  xlab("Concentration (µM)")+
  scale_x_continuous(breaks = c(0,6.25,12.5,25,50,100),labels = c("0","6.25","12.5","25","50","100"),expand = c(0,0.5))+
  scale_fill_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V"), values = cols)+
  scale_color_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V"), values = cols)+
  theme(axis.title.x = element_text(size=12, face = "bold",margin = unit(c(0 ,0,0,0),units = "in")),
        axis.text.x = element_text(size=12,face = "bold",angle = 90,vjust = 0.5,hjust=0.5,margin = unit(c(0,0,0,0),units = "in")),
        axis.text.y = element_text(size=12, face = "bold",margin = unit(c(0 ,0,0,0),units = "in")),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size = 0.5),
        axis.title.y=element_text(size=12,face = "bold"),
        legend.key.size = unit(1,"in"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        plot.margin = unit(c(0, 0, 0, 0), "in"),
        legend.position = "None")

#Figure 1 D

load("S5.RData")

fbz_stats <- aov(phenotype ~ strain, data = regressed_fbzH)%>%
  rstatix::tukey_hsd()%>%
  dplyr::mutate(p.adj.signif = case_when(p.adj > 0.05 ~ "ns",
                                         (p.adj < 0.05 & p.adj > 0.01) ~"*",
                                         (p.adj < 0.01 & p.adj > 0.001) ~"**",
                                         (p.adj < 0.001 & p.adj > 0.0001) ~"***",
                                         (p.adj < 0.0001 & p.adj > 0.00001) ~"****",
                                         (p.adj < 0.00001) ~ "*****"))%>%
  dplyr::filter(grepl("N2",group1) | grepl("N2",group2))


regressedfbz_medianEXTold <-regressed_fbzH%>%
  dplyr::filter(!is.na(strain))%>%
  dplyr::filter(trait == traita)%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "882","919","920","1081","1082","1075","1076","1327","1328", "1325","1326")))%>%
  ggplot()+
  aes(x=fancy_strain, y = phenotype)+
  geom_jitter(width = 0.1,size=0.1)+
  geom_boxplot(aes(fill=strain, alpha=0.1), outlier.shape = NA)+
  theme_cowplot(8)+
  ylab(glue('Normalized OD'))+
  scale_y_continuous(limits = c(-400,250),sec.axis = dup_axis(name = "Fenbendazole"))+
  stat_pvalue_manual(fbz_stats, label = "p.adj.signif",xmax="group1", y.position = c(225),remove.bracket = TRUE,size = 4)+
  scale_x_discrete(labels=c("N2" = "N2", "882" = "Deletion","920" = "F200Y","1082"="E198A","1327"="E198V","1325"="E198L","1076"="F167Y"))+
  scale_fill_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V"), values = cols)+
  scale_color_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V"), values = cols)+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(size=12,face = "bold",angle = 45,hjust = 1,vjust = 1,margin = unit(c(0,0,0,0),units = "in")),
        axis.text.y = element_text(size=12, face = "bold",margin = unit(c(0 ,0,0,0),units = "in")),
        axis.title.y=element_text(size = 12,face = "bold",margin = unit(c(0,0.075,0,0),units = "in")),
        legend.key.size = unit(0.4,"in"),
        axis.text.y.right = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y.right = element_text(margin = unit(c(0,0,0,0),units = "in"), angle=-90, hjust = 0.5,vjust = 0.5),
        axis.ticks.y.right = element_blank(),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0), "in"),
        legend.position = "None")
# Plot and save full figure 
con_alleles <- cowplot::plot_grid(abzmedian.EXT_subtractedold,regressedabz_medianEXTold,FBZmedian.EXT_subtractedold,regressedfbz_medianEXTold,ncol=2,nrow=2,labels = c("A","B","C","D"),align = "vh",axis = "lrbt",label_size = 12,label_fontfamily = "Arial",rel_widths = c(1,1,1,1),rel_heights = c(1,1,1,1))
ggsave("~/Desktop/ben1_2020_CMD/manuscript/figure_1.jpeg",plot = con_alleles,device = "jpeg",width = 7.5,height = 6,units = "in")

##Figure 2

load("S1.RData")
abz_outpruned <- subtracted_dose

abz_outprunedmut <- abz_outpruned%>%
  dplyr::filter(condition != "200")%>%
  dplyr::filter(trait == traita)%>%
  dplyr::filter(strain %in% para_new)%>%
  dplyr::mutate(numcon = case_when(condition == "DMSO" ~ 0,
                                   condition == "6_25" ~ 6.25,
                                   condition == "12_5" ~ 12.5,
                                   condition == "25" ~ 25,
                                   condition == "50" ~ 50,
                                   condition == "100" ~ 100,
                                   condition == "200" ~ 200))%>%
  dplyr::group_by(strain, condition)%>%
  dplyr::mutate(meancondition = mean(sub_cntrl))%>%
  dplyr::mutate(SND = sd(sub_cntrl))

abzmedian.EXT_subtractednew <- abz_outprunedmut%>%
  dplyr::filter(!is.na(strain))%>%
  dplyr::filter(strain %in% c("N2","919","920","882","1081","1082", "1075","1076","1325","1326","1327","1328"))%>%
  dplyr::filter(condition != "200")%>%
  dplyr::filter(numcon %in% c(0,6.25,12.5,25,50,100,200))%>%
  dplyr::filter(trait == traita)%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "882","919","920","1081","1082","1327","1328", "1325","1326","1075","1076")))%>%
  dplyr::mutate(con_fix = factor(condition, levels = c("DMSO","6_25","12_5","25","50","100","200")))%>%
  ggplot()+
  aes(x=numcon, y = meancondition)+
  geom_errorbar(aes(color = fancy_strain, ymin=meancondition - SND, ymax=meancondition + SND), width=2, size = 1)+
  geom_line(aes(x = numcon, y = meancondition, colour= fancy_strain), size = 1)+
  theme_cowplot(12)+
  scale_y_continuous(limits = c(-600,150))+
  ylab("Difference in OD")+
  xlab("Concentration (µM)")+
  scale_x_continuous(breaks = c(0,6.25,12.5,25,50,100),labels = c("0","6.25","12.5","25","50","100"),expand = c(0,0.5))+
  scale_fill_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V"), values = cols)+
  scale_color_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V"), values = cols)+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(size=12,face = "bold",angle = 90,vjust = 0.5,hjust=0.5,margin = unit(c(0,0,0,0),units = "in")),
        axis.text.y = element_text(size=12, face = "bold"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size = 0.5),
        axis.title.y=element_text(size=12,face = "bold"),
        legend.key.size = unit(1,"in"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        plot.margin = unit(c(0, 0, 0, 0), "in"),
        legend.position = "None")
#Figure 2 B

load("S3.RData")
abz_stats <- aov(phenotype ~ strain, data = regressedabz_abzHN)%>%
  rstatix::tukey_hsd()%>%
  dplyr::mutate(p.adj.signif = case_when(p.adj > 0.05 ~ "ns",
                                         (p.adj < 0.05 & p.adj > 0.01) ~"*",
                                         (p.adj < 0.01 & p.adj > 0.001) ~"**",
                                         (p.adj < 0.001 & p.adj > 0.0001) ~"***",
                                         (p.adj < 0.0001 & p.adj > 0.00001) ~"****",
                                         (p.adj < 0.00001) ~ "*****"))%>%
  dplyr::filter(grepl("N2",group1) | grepl("N2",group2))


regressedabz_medianEXTnew <-regressedabz_abzHN%>%
  dplyr::filter(!is.na(strain))%>%
  dplyr::filter(trait == traita)%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "882","919","920","1081","1082","1075","1076","1327","1328", "1325","1326")))%>%
  ggplot()+
  aes(x=fancy_strain, y = phenotype)+
  geom_jitter(width = 0.1,size=0.1)+
  geom_boxplot(aes(fill=strain, alpha = 0.1), outlier.shape = NA)+
  theme_cowplot(8)+
  scale_y_continuous(limits = c(-400, 250), sec.axis = dup_axis(name = "Albendazole"))+
  stat_pvalue_manual(abz_stats, label = "p.adj.signif",xmax="group1", y.position = c(220),remove.bracket = TRUE,size = 4)+
  ylab(glue("Normalized OD"))+
  scale_x_discrete(labels=c("N2" = "N2", "882" = "Deletion","919" = "F200Y","1082"="E198A","1327"="E198V","1325"="E198L","1076"="F167Y"))+
  scale_fill_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V"), values = cols)+
  scale_color_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V"), values = cols)+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=12, face = "bold",margin = unit(c(0 ,0,0,0),units = "in")),
        axis.title.y=element_text(size = 12,face = "bold",margin = unit(c(0,0.075,0,0),units = "in")),
        legend.key.size = unit(0.4,"in"),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(margin = unit(c(0,0,0,0),units = "in")),
        axis.ticks.y.right = element_blank(),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        plot.margin = unit(c(0.1, 0.1, 0, 0), "in"),
        legend.position = "None")

#Figure 2 C


load("S4.RData")
fbz_outpruned <- subtracted_dose_fbz

fbz_outprunedmut <- fbz_outpruned%>%
  dplyr::filter(condition != "200")%>%
  dplyr::filter(trait == traita)%>%
  dplyr::mutate(numcon = case_when(condition == "DMSO" ~ 0,
                                   condition == "6_25" ~ 6.25,
                                   condition == "12_5" ~ 12.5,
                                   condition == "25" ~ 25,
                                   condition == "50" ~ 50,
                                   condition == "100" ~ 100))%>%
  dplyr::group_by(strain, condition)%>%
  dplyr::mutate(meancondition = mean(sub_cntrl))%>%
  dplyr::mutate(SND = sd(sub_cntrl))

FBZmedian.EXT_subtractednew <-fbz_outprunedmut%>%
  dplyr::filter(!is.na(strain))%>%
  dplyr::filter(strain %in% c("N2","882","1325", "1327"))%>%
  dplyr::filter(condition != "200")%>%
  dplyr::filter(numcon %in% c(0,6.25,12.5,25,50,100))%>%
  dplyr::filter(trait == traita)%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "882","919","920","1081","1082","1327","1328", "1325","1326","1075","1076")))%>%
  dplyr::mutate(con_fix = factor(condition, levels = c("DMSO","6_25","12_5","25","50","100")))%>%
  ggplot()+
  aes(x=numcon, y = meancondition)+
  geom_errorbar(aes(color = fancy_strain, ymin=meancondition - SND, ymax=meancondition + SND), width=2, size = 1)+
  geom_line(aes(x = numcon, y = meancondition, colour= fancy_strain), size = 1)+
  theme_cowplot(12)+
  scale_y_continuous(limits = c(-750,100))+
  ylab("Difference in OD")+
  xlab("Concentration (µM)")+
  scale_x_continuous(breaks = c(0,6.25,12.5,25,50,100),labels = c("0","6.25","12.5","25","50","100"),expand = c(0,0.5))+
  scale_fill_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V"), values = cols)+
  scale_color_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V"), values = cols)+
  theme(axis.title.x = element_text(size=12, face = "bold",margin = unit(c(0 ,0,0,0),units = "in")),
        axis.text.x = element_text(size=12,face = "bold",angle = 90,vjust = 0.5,hjust=0.5,margin = unit(c(0,0,0,0),units = "in")),
        axis.text.y = element_text(size=12, face = "bold"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size = 0.5),
        axis.title.y=element_text(size=12,face = "bold"),
        legend.key.size = unit(1,"in"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        plot.margin = unit(c(0, 0, 0, 0), "in"),
        legend.position = "None")

#Figure 2 D

load("S6.RData")

fbz_stats <- aov(phenotype ~ strain, data = regressed_fbzHN)%>%
  rstatix::tukey_hsd()%>%
  dplyr::mutate(p.adj.signif = case_when(p.adj > 0.05 ~ "ns",
                                         (p.adj < 0.05 & p.adj > 0.01) ~"*",
                                         (p.adj < 0.01 & p.adj > 0.001) ~"**",
                                         (p.adj < 0.001 & p.adj > 0.0001) ~"***",
                                         (p.adj < 0.0001 & p.adj > 0.00001) ~"****",
                                         (p.adj < 0.00001) ~ "*****"))%>%
  dplyr::filter(grepl("N2",group1) | grepl("N2",group2))


regressedfbz_medianEXTnew <-regressed_fbzHN%>%
  dplyr::filter(!is.na(strain))%>%
  dplyr::filter(trait == traita)%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "882","919","920","1081","1082","1075","1076","1327","1328", "1325","1326")))%>%
  ggplot()+
  aes(x=fancy_strain, y = phenotype)+
  geom_jitter(width = 0.1,size=0.1)+
  geom_boxplot(aes(fill=strain, alpha=0.1), outlier.shape = NA)+
  theme_cowplot(8)+
  ylab(glue('Normalized OD'))+
  scale_y_continuous(limits = c(-400,250),sec.axis = dup_axis(name = "Fenbendazole"))+
  stat_pvalue_manual(fbz_stats, label = "p.adj.signif",xmax="group1", y.position = c(225),remove.bracket = TRUE,size = 4)+
  scale_x_discrete(labels=c("N2" = "N2", "882" = "Deletion","920" = "F200Y","1082"="E198A","1327"="E198V","1325"="E198L","1076"="F167Y"))+
  scale_fill_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V"), values = cols)+
  scale_color_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V"), values = cols)+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(size=12,face = "bold",angle = 45,hjust = 1,vjust = 1,margin = unit(c(0,0,0,0),units = "in")),
        axis.text.y = element_text(size=12, face = "bold",margin = unit(c(0 ,0,0,0),units = "in")),
        axis.title.y=element_text(size = 12,face = "bold",margin = unit(c(0,0.075,0,0),units = "in")),
        legend.key.size = unit(0.4,"in"),
        axis.text.y.right = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y.right = element_text(margin = unit(c(0,0,0,0),units = "in"), angle=-90, hjust = 0.5,vjust = 0.5),
        axis.ticks.y.right = element_blank(),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0), "in"),
        legend.position = "None")
## Plot and save

new_alleles <- cowplot::plot_grid(abzmedian.EXT_subtractednew,regressedabz_medianEXTnew,FBZmedian.EXT_subtractednew,regressedfbz_medianEXTnew,ncol=2,nrow=2,labels = c("A","B","C","D"),align = "hv",label_size = 12,label_fontfamily = "Arial",rel_widths = c(1,1,1,1),rel_heights = c(1,1,1,1))
ggsave("~/Desktop/ben1_2020_CMD/manuscript/figure_2.jpeg",plot = new_alleles,device = "jpeg",width = 7.5,height = 6,units = "in")


## Figure 3

#Figure 3 A
comp_data <- readr::read_tsv("S10.tsv")
comp_data$Allele_Freq <- as.numeric(comp_data$Allele_Freq)
comp_data$Mean <- as.numeric(comp_data$Mean)
comp_data <- mutate(comp_data, Strain = case_when(Strain == "del" ~ "Deletion",
                    TRUE ~ as.character(Strain)))
cols <- c("N2" = "orange", "Deletion" = "grey", "F200Y" = "Red", "E198V"= "Purple","E198L"="yellow","E198A"="blue","F167Y"="green")
x <- comp_data %>%
  dplyr::select(Strain,Allele_Freq,Condition,Generation,Rep,Mean)%>%
  na.omit()%>%
  dplyr::group_by(Condition, Strain, Generation)%>%
  dplyr::mutate(Mean = mean(1 - Allele_Freq)) %>%
  dplyr::mutate(SND = sd(1-Allele_Freq)) %>%
  dplyr::ungroup()%>%
  dplyr::mutate(Generation = case_when(Generation == 1 ~ 1,
                                       Generation == 2 ~ 3,
                                       Generation == 3 ~ 5,
                                       Generation == 4 ~ 7))%>%
  dplyr::mutate(Condition = case_when(Condition == "A" ~ "Albendazole",
                                      Condition == "D" ~ "DMSO"))
x$Condition_f = factor(x$Condition, levels = c("DMSO","Albendazole"))

DMSO_comp <- x %>%
  dplyr::filter(Condition_f == "DMSO")%>%
  ggplot()+
  aes(x = Generation, y = Mean , fill = Condition_f, color = Strain) +
  geom_errorbar(aes(ymin=Mean - SND, ymax=Mean + SND), width=0.25,size=0.5)+
  geom_line(aes(linetype = Condition_f),size=1)+
  scale_color_manual(labels = c("N2"="N2","del" = "Deletion","F200Y"="F200Y","E198V"="E198V","E198L"="E198L","E198A" = "E198A","F167Y"="F167Y"), values = cols)+
  theme_cowplot(8)+
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),limits = c(0,1.0))+
  scale_x_continuous(breaks = c(1,3,5,7))+
  ylab("Relative allele frequency")+
  xlab("Generation")+
  theme(axis.title.x=element_text(size = 12,face = "bold"),
        axis.text.x = element_text(size = 12,face = "bold"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text.y=element_text(size=12, face = "bold",margin = unit(c(0 ,0,0,0),units = "in")),
        axis.title.y=element_text(size = 12, face = "bold",margin = unit(c(0,0.075,0,0),units = "in")),
        legend.key.size = unit(0.4,"in"),
        plot.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.text = element_text(face = "bold"),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "in"),
        legend.position = "None")
#Figure 3 B
fitness <- readr::read_tsv("S11.txt")
fitness <- mutate(fitness, Strain = case_when(Strain == "del" ~ "Deletion",
                                                  TRUE ~ as.character(Strain)))
Dfitness <- dplyr::filter(fitness,Condition == "D")
dmso_stats <- aov(Fitness ~ Strain, data = Dfitness)%>%
  rstatix::tukey_hsd()%>%
  dplyr::mutate(p.adj.signif = case_when(p.adj > 0.05 ~ "ns",
                                         (p.adj < 0.05 & p.adj > 0.01) ~"*",
                                         (p.adj < 0.01 & p.adj > 0.001) ~"**",
                                         (p.adj < 0.001 & p.adj > 0.0001) ~"***",
                                         (p.adj < 0.0001 & p.adj > 0.00001) ~"****",
                                         (p.adj < 0.00001) ~ "*****"))%>%
  dplyr::filter(grepl("N2",group1) | grepl("N2",group2))
DFIT <- fitness %>%
  dplyr::filter(Condition == "D")%>%
  dplyr::mutate(fancy_strain=factor(Strain, levels = c("N2", "Deletion","F200Y","E198A","F167Y","E198V","E198L")))%>%
  ggplot()+
  aes(x=fancy_strain, y=Fitness)+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_jitter(width = 0.1)+
  geom_boxplot(aes(fill = Strain, alpha=0.1),outlier.shape = NA)+
  theme_cowplot(8)+
  stat_pvalue_manual(dmso_stats, label = "p.adj.signif",xmax="group1", y.position = c(0.15),remove.bracket = TRUE)+
  ylab("Competitive Fitness")+
  scale_x_discrete(labels = c("N2"="N2","del" = "Deletion","F200Y"="F200Y","E198A" ="E198A","F167Y"="F167Y","E198L"="E198L","E198V"="E198V"))+
  scale_color_manual(labels = c("N2"="N2","del" = "Deletion","F200Y"="F200Y","E198V"="E198V"), values = cols)+
  scale_fill_manual(labels = c("N2"="N2","del" = "Deletion","F200Y"="F200Y","E198V"="E198V"), values = cols)+
  theme(legend.position="none")+
  xlab("Strain")+
  scale_y_continuous(sec.axis = dup_axis(name = "DMSO"))+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(size = 12,face = "bold", angle = 45,hjust = 1,vjust = 1,margin = unit(c(0,0,0,0),units = "in")),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y=element_text(size=12, face = "bold",margin = unit(c(0 ,0,0,0),units = "in")),
        axis.title.y=element_text(size = 12, face = "bold",margin = unit(c(0,0.075,0,0),units = "in")),
        legend.key.size = unit(0.4,"in"),
        legend.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.text = element_text(face = "bold"),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "in"),
        legend.position = "None")

#Figure 3 C

ABZ_comp <- x %>%
  dplyr::filter(Condition_f == "Albendazole")%>%
  ggplot()+
  aes(x = Generation, y = Mean , fill = Condition_f, color = Strain) +
  geom_errorbar(aes(ymin=Mean - SND, ymax=Mean + SND), width=0.25,size=0.5)+
  geom_line(aes(linetype = Condition_f), size=1)+ 
  scale_color_manual(labels = c("N2"="N2","del" = "Deletion","F200Y"="F200Y","E198V"="E198V","E198L"="E198L","E198A" = "E198A","F167Y"="F167Y"), values = cols)+
  theme_cowplot(8)+
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),limits = c(0,1.0))+
  scale_x_continuous(breaks = c(1,3,5,7))+
  ylab("Relative allele frequency")+
  xlab("Generation")+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(size = 12,face = "bold"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y=element_text(size=12, face = "bold",margin = unit(c(0 ,0,0,0),units = "in")),
        axis.title.y=element_text(size = 12, face = "bold",margin = unit(c(0,0.075,0,0),units = "in")),
        legend.key.size = unit(0.4,"in"),
        #axis.ticks.x = element_blank(),
        plot.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.text = element_text(face = "bold"),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "in"),
        legend.position = "None")
#Figure 3 D

Afitness <- dplyr::filter(fitness,Condition == "A")
A_stats <- aov(Fitness ~ Strain, data = Afitness)%>%
  rstatix::tukey_hsd()%>%
  dplyr::mutate(p.adj.signif = case_when(p.adj > 0.05 ~ "ns",
                                         (p.adj < 0.05 & p.adj > 0.01) ~"*",
                                         (p.adj < 0.01 & p.adj > 0.001) ~"**",
                                         (p.adj < 0.001 & p.adj > 0.0001) ~"***",
                                         (p.adj < 0.0001 & p.adj > 0.00001) ~"****",
                                         (p.adj < 0.00001) ~ "*****"))%>%
  dplyr::filter(grepl("N2",group1) | grepl("N2",group2))
A_F200YE198Vstats <- aov(Fitness ~ Strain, data = Afitness)%>%
  rstatix::tukey_hsd()%>%
  dplyr::mutate(p.adj.signif = case_when(p.adj > 0.05 ~ "ns",
                                         (p.adj < 0.05 & p.adj > 0.01) ~"*",
                                         (p.adj < 0.01 & p.adj > 0.001) ~"**",
                                         (p.adj < 0.001 & p.adj > 0.0001) ~"***",
                                         (p.adj < 0.0001 & p.adj > 0.00001) ~"****",
                                         (p.adj < 0.00001) ~ "*****"))%>%
  dplyr::filter(grepl("E198V",group1) , grepl("F200Y",group2))
AFIT <- fitness %>%
  dplyr::filter(Condition == "A")%>%
  dplyr::mutate(fancy_strain=factor(Strain, levels = c("N2", "Deletion","F200Y","E198A","F167Y","E198V","E198L")))%>%
  ggplot()+
  aes(x=fancy_strain, y=Fitness)+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_jitter(width = 0.1)+
  geom_boxplot(aes(fill = Strain, alpha=0.1),outlier.shape = NA)+
  theme_cowplot(8)+
  stat_pvalue_manual(A_stats, label = "p.adj.signif",xmax="group1", y.position = c(0.75),remove.bracket = TRUE)+
  stat_pvalue_manual(A_F200YE198Vstats, label = "p.adj.signif",y.position = c(0.70), remove.bracket = FALSE)+
  scale_x_discrete(labels = c("N2"="N2","del" = "Deletion","F200Y"="F200Y","E198A" ="E198A","F167Y"="F167Y","E198L"="E198L","E198V"="E198V"))+
  scale_color_manual(labels = c("N2"="N2","del" = "Deletion","F200Y"="F200Y","E198V"="E198V"), values = cols)+
  scale_fill_manual(labels = c("N2"="N2","del" = "Deletion","F200Y"="F200Y","E198V"="E198V"), values = cols)+
  theme(legend.position="none")+
  scale_y_continuous(limits = c(-0.2,0.8),sec.axis = dup_axis(name = "Albendazole"))+
  ylab("Competitive Fitness")+
  xlab("Strain")+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(size = 12,face = "bold", angle = 45,hjust = 1,vjust = 1,margin = unit(c(0,0,0,0),units = "in")),
        axis.text.y.right = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y=element_text(size=12, face = "bold",margin = unit(c(0 ,0,0,0),units = "in")),
        axis.title.y=element_text(size = 12, face = "bold",margin = unit(c(0,0.075,0,0),units = "in")),
        legend.key.size = unit(0.4,"in"),
        legend.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.text = element_text(face = "bold"),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "in"),
        legend.position = "None")

##Plot and Save
competitionplot <- cowplot::plot_grid(ABZ_comp,AFIT,DMSO_comp,DFIT,nrow = 2,ncol = 2,labels = c("A","B","C","D"),align = "vh",axis = "lrbt",label_size = 12,label_fontfamily = "Arial",rel_widths = c(1,1,1,1),rel_heights = c(1,1,1,1))
ggsave("~/Desktop/ben1_2020_CMD/manuscript/figure_3.jpeg",plot = competitionplot,device = "jpeg",width = 7.5,height = 6,units = "in")

### SUPPLEMENTAL FIGURES

## Supplemental Figure 2
load("S1.RData")
abz_outpruned <- subtracted_dose

abz_outprunedmut <- abz_outpruned%>%
  dplyr::filter(condition != "200")%>%
  dplyr::filter(trait == traita)%>%
  dplyr::mutate(numcon = case_when(condition == "DMSO" ~ 0,
                                   condition == "6_25" ~ 6.25,
                                   condition == "12_5" ~ 12.5,
                                   condition == "25" ~ 25,
                                   condition == "50" ~ 50,
                                   condition == "100" ~ 100,
                                   condition == "200" ~ 200))%>%
  dplyr::group_by(strain, condition)%>%
  dplyr::mutate(meancondition = mean(sub_cntrl))%>%
  dplyr::mutate(SND = sd(sub_cntrl))

load("S4.RData")
fbz_outpruned <- subtracted_dose_fbz

fbz_outprunedmut <- fbz_outpruned%>%
  dplyr::filter(condition != "200")%>%
  dplyr::filter(trait == traita)%>%
  dplyr::mutate(numcon = case_when(condition == "DMSO" ~ 0,
                                   condition == "6_25" ~ 6.25,
                                   condition == "12_5" ~ 12.5,
                                   condition == "25" ~ 25,
                                   condition == "50" ~ 50,
                                   condition == "100" ~ 100))%>%
  dplyr::group_by(strain, condition)%>%
  dplyr::mutate(meancondition = mean(sub_cntrl))%>%
  dplyr::mutate(SND = sd(sub_cntrl))


fbz_strainsplit <- fbz_outprunedmut%>%
  dplyr::filter(!is.na(strain))%>%
  dplyr::filter(strain %in% plate_1strains)%>%
  dplyr::filter(condition != "200")%>%
  dplyr::filter(trait == traita)%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "882","919","920","1081","1082","1327","1328", "1325","1326","1075","1076")))%>%
  dplyr::mutate(fixedname = case_when(fancy_strain %in% c("N2") ~ "N2",
                                      fancy_strain %in% c("882") ~ "Deletion",
                                      fancy_strain %in% c("919","920") ~ "F200Y",
                                      fancy_strain %in% c("1081", "1082") ~ "E198A",
                                      fancy_strain %in% c("1327","1328") ~ "E198V",
                                      fancy_strain %in% c("1325","1326") ~ "E198L",
                                      fancy_strain %in% c("1075","1076") ~ "F167Y"))%>%
  dplyr::mutate(fixedname = factor(fixedname, levels = c("N2", 'Deletion',"F200Y","E198A",'E198V',"E198L",'F167Y')))%>%
  dplyr::mutate(con_fix = factor(condition, levels = c("DMSO","6_25","12_5","25","50","100")))%>%
  dplyr::filter(fixedname != "N2")%>%
  dplyr::filter(fixedname != "Deletion")%>%
  ggplot()+
  geom_errorbar(aes(color = fancy_strain, ymin=meancondition - SND, ymax=meancondition + SND), width=3, size = 0.25)+
  aes(x=numcon, y = phenotype)+
  geom_line(aes(x = numcon, y = meancondition, colour = fancy_strain), size = 0.5)+
  theme_cowplot(12)+
  xlab("Fenbendazole concentration (µM)")+
  ylab("Normalized OD")+
  ylim(-850,150)+
  scale_x_continuous(breaks = c(0,6.25,12.5,25,50,100),expand = c(0,0.5),labels = c("0","6.25","12.5","25","50","100"))+
  scale_fill_manual(name = "fancy_strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V"), values = c("N2" = "orange","882"="grey","919"="red","920"="red","1325"="yellow","1326"="yellow","1327"="purple","1328"="purple","1082" = "blue","1081" = "blue", "1076" = "green","1075" = "green"))+
  scale_color_manual(name = "fancy_strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V"), values = c("N2" = "orange","882"="grey","919"="red","920"="red","1325"="yellow","1326"="yellow","1327"="purple","1328"="purple","1082" = "blue","1081" = "blue", "1076" = "green","1075" = "green"))+
  facet_grid(rows = vars(fixedname))+
  theme(axis.title.x=element_text(face = "bold", size = 12, margin = unit(c(0.1,0,0,0), units = "in"),vjust=1),
        axis.text.x = element_text(size = 12,face = "bold", margin = unit(c(0.05,0,0,0), units = "in"), hjust = 1, vjust = 0.5,angle=90),
        axis.text.y = element_text(size=12, face = "bold",margin = unit(c(0 ,0,0,0),units = "in")),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.ticks.x = element_line(size = 0.25),
        axis.ticks.y = element_line(size = 0.5),
        axis.title.y=element_blank(),
        strip.background = element_blank(),
        axis.title.y.right=element_text(size=12,face = "bold", margin = unit(c(0,0,0,0.01),units = "in")),
        legend.key.size = unit(0.4,"in"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        plot.margin = unit(c(t=0.2, 0, 0, 0), "in"), 
        legend.position = "None")

abz_strainsplit <- abz_outprunedmut%>%
  dplyr::filter(!is.na(strain))%>%
  dplyr::filter(strain %in% plate_1strains)%>%
  dplyr::filter(condition != "200")%>%
  dplyr::filter(trait == traita)%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "882","919","920","1081","1082","1327","1328", "1325","1326","1075","1076")))%>%
  dplyr::mutate(fixedname = case_when(fancy_strain %in% c("N2") ~ "N2",
                                      fancy_strain %in% c("882") ~ "Deletion",
                                      fancy_strain %in% c("919","920") ~ "F200Y",
                                      fancy_strain %in% c("1081", "1082") ~ "E198A",
                                      fancy_strain %in% c("1327","1328") ~ "E198V",
                                      fancy_strain %in% c("1325","1326") ~ "E198L",
                                      fancy_strain %in% c("1075","1076") ~ "F167Y"))%>%
  dplyr::mutate(fixedname = factor(fixedname, levels = c("N2", 'Deletion',"F200Y","E198A",'E198V',"E198L",'F167Y')))%>%
  dplyr::mutate(con_fix = factor(condition, levels = c("DMSO","6_25","12_5","25","50","100")))%>%
  dplyr::filter(fixedname != "N2")%>%
  dplyr::filter(fixedname != "Deletion")%>%
  ggplot()+
  geom_errorbar(aes(color = fancy_strain, ymin=meancondition - SND, ymax=meancondition + SND), width=3, size = 0.25)+
  aes(x=numcon, y = phenotype)+
  geom_line(aes(x = numcon, y = meancondition, colour = fancy_strain), size = 0.5)+
  theme_cowplot(12)+
  xlab("Albendazole concentration (µM)")+
  ylab("Difference in OD")+
  ylim(-500,150)+
  scale_x_continuous(breaks = c(0,6.25,12.5,25,50,100),expand = c(0,0.5),labels = c("0","6.25","12.5","25","50","100"))+
  scale_fill_manual(name = "fancy_strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V"), values = c("N2" = "orange","882"="grey","919"="red","920"="red","1325"="yellow","1326"="yellow","1327"="purple","1328"="purple","1082" = "blue","1081" = "blue", "1076" = "green","1075" = "green"))+
  scale_color_manual(name = "fancy_strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V"), values = c("N2" = "orange","882"="grey","919"="red","920"="red","1325"="yellow","1326"="yellow","1327"="purple","1328"="purple","1082" = "blue","1081" = "blue", "1076" = "green","1075" = "green"))+
  facet_grid(rows = vars(fixedname))+
  theme(axis.title.x=element_text(face = "bold", size = 12, margin = unit(c(0.1,0,0,0), units = "in"),vjust=1),
        axis.text.x = element_text(size = 12,face = "bold", margin = unit(c(0.05,0,0,0), units = "in"), hjust = 1, vjust = 0.5,angle=90),
        axis.text.y = element_text(size=12, face = "bold",margin = unit(c(0 ,0,0,0),units = "in")),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.ticks.x = element_line(size = 0.25),
        axis.ticks.y = element_line(size = 0.5),
        axis.title.y=element_text(size=12,face = "bold"),
        axis.title.y.right=element_text(size=12,face = "bold", margin = unit(c(0,0,0,0.01),units = "in")),
        legend.key.size = unit(0.4,"in"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        legend.position = "None")
gridallele <- cowplot::plot_grid(abz_strainsplit,fbz_strainsplit, nrow = 1, ncol = 2,labels = c("A","B"),align = "vh",axis = "lrbt",label_size = 12,label_fontfamily = "Arial",rel_widths = c(1,1,1,1),rel_heights = c(1,1,1,1))
ggsave(filename = "~/Desktop/ben1_2020_CMD/manuscript/Supplemental_figure_3.jpeg", plot = gridallele, device = "jpeg",units = "in",width = 7,height = 5)

##Supplemental Figure 3 

load("S7.RData")
allele_Deletion <- rio::import("S12.txt")
allele_Deletion <- allele_Deletion %>%
  dplyr::mutate(., group1 = as.character(group1))%>%
  dplyr::mutate(., group2 = as.character(group2))
bothfbz_stats <- aov(phenotype ~ strain, data = regressed_fbzHA)%>%
  rstatix::tukey_hsd()%>%
  dplyr::mutate(p.adj.signif = case_when(p.adj > 0.05 ~ "ns",
                                         (p.adj < 0.05 & p.adj > 0.01) ~"*",
                                         (p.adj < 0.01 & p.adj > 0.001) ~"**",
                                         (p.adj < 0.001 & p.adj > 0.0001) ~"***",
                                         (p.adj < 0.0001 & p.adj > 0.00001) ~"****",
                                         (p.adj < 0.00001) ~ "*****"))%>%
  dplyr::filter(grepl("N2",group1) | grepl("N2",group2))

regressedfbz_medianEXT_bothalleles <-regressed_fbzHA%>%
  dplyr::filter(!is.na(strain))%>%
  dplyr::filter(trait == traita)%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "882","920","1081","1082","1075","1076","1327","1328", "1325","1326")))%>%
  ggplot()+
  aes(x=fancy_strain, y = phenotype)+
  geom_jitter(width = 0.1,size=0.1)+
  geom_boxplot(aes(fill=strain, alpha=0.1), outlier.shape = NA)+
  stat_pvalue_manual(bothfbz_stats, label = "p.adj.signif",xmax="group1", y.position = c(170),remove.bracket = TRUE)+
  stat_pvalue_manual(allele_Deletion, label = "strain", xmax = "group2", y.position = c(220), remove.bracket = TRUE)+
  theme_cowplot(8)+
  ylab(glue('Normalized OD'))+
  scale_x_discrete(labels=c("N2" = "N2", "882" = expression(paste(italic("ean64"))),"920" = expression(italic("ean101")),"1081"=expression(italic("ean149")),"1082"=expression(italic("ean150")),"1327"=expression(italic("ean203")),"1328"=expression(italic("ean204")),"1325"=expression(italic("ean201")),"1326"=expression(italic("ean202")),"1075"=expression(italic("ean143")),"1076" = expression(italic("ean144"))))+
  scale_fill_manual(name = "fancy_strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V"), values = c("N2" = "orange","882"="grey","919"="red","920"="red","1325"="yellow","1326"="yellow","1327"="purple","1328"="purple","1082" = "blue","1081" = "blue", "1076" = "green","1075" = "green"))+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(size=12,face = "bold",angle = 45,hjust = 1,vjust = 1,margin = unit(c(0,0,0,0),units = "in")),
        axis.text.y = element_text(size=12, face = "bold"),
        axis.title.y=element_text(size = 12,face = "bold"),
        legend.key.size = unit(0.4,"in"),
        axis.text.y.right = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y.right = element_text(margin = unit(c(0,0,0,0),units = "in"), hjust = 0.1),
        axis.ticks.y.right = element_blank(),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        legend.position = "None")

ggsave("~/Desktop/ben1_2020_CMD/manuscript/Supplemental_figure_4.jpeg",plot = regressedfbz_medianEXT_bothalleles,device = "jpeg",width = 7.5,height = 4,units = "in")

##Supplemental Figure 4

load("S8.RData")
DMSOabz <- biopruned_abzH %>%
  dplyr::filter(strain %in% c("N2","882","919", "1076", "1082","1325", "1327"))%>%
  dplyr::filter(condition == "DMSO")

dmsoabz_stats <- aov(median.EXT ~ strain, data = DMSOabz)%>%
  rstatix::tukey_hsd()%>%
  dplyr::mutate(p.adj.signif = case_when(p.adj > 0.05 ~ "ns",
                                         (p.adj < 0.05 & p.adj > 0.01) ~"*",
                                         (p.adj < 0.01 & p.adj > 0.001) ~"**",
                                         (p.adj < 0.001 & p.adj > 0.0001) ~"***",
                                         (p.adj < 0.0001 & p.adj > 0.00001) ~"****",
                                         (p.adj < 0.00001) ~ "*****"))%>%
  dplyr::filter(grepl("N2",group1) | grepl("N2",group2))

DMSOabz_medianEXT <-DMSOabz%>%
  dplyr::filter(!is.na(strain))%>%
  dplyr::filter(strain %in% parafull)%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "882","919","920","1081","1082","1075","1076","1327","1328", "1325","1326")))%>%
  ggplot()+
  aes(x=fancy_strain, y = median.EXT)+
  geom_jitter(width = 0.1,size=0.1)+
  geom_boxplot(aes(fill=strain,alpha=0.1), outlier.shape = NA)+
  theme_cowplot(8)+
  ylab(glue('OD'))+
  stat_pvalue_manual(dmsoabz_stats, label = "p.adj.signif",xmax="group1", y.position = c(1550),remove.bracket = TRUE)+
  scale_x_discrete(labels=c("N2" = "N2", "882" = "Deletion","919" = "F200Y","1082"="E198A","1327"="E198V","1325"="E198L","1076"="F167Y"))+
  scale_fill_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","F200Y"="F200Y", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V"), values = c("N2" = "orange","882"="grey","919"="red","920"="red","1325"="yellow","1326"="yellow","1327"="purple","1328"="purple","1082" = "blue","1081" = "blue", "1076" = "green","1075" = "green"))+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(size=12, face = "bold",margin = unit(c(0 ,0,0,0),units = "in")),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=12, face = "bold",margin = unit(c(0 ,0,0,0),units = "in")),
        axis.title.y=element_text(size = 12,face = "bold",margin = unit(c(0,0.075,0,0),units = "in")),
        legend.key.size = unit(0.4,"in"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        legend.position = "None")

ggsave("~/Desktop/ben1_2020_CMD/manuscript/Supplemental_figure_2.jpeg",plot = DMSOabz_medianEXT,device = "jpeg",width = 7.5,height = 4,units = "in")


##Supplemental Figure 5

load("S4.RData")
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
save(collectEC, file = "S9.RData")
allele_F <- rio::import("S12.txt")
allele_F <- allele_F %>%
  dplyr::mutate(group1 = case_when(group1 == "920" ~ "919",
                TRUE ~ as.character(group1)))%>%
  dplyr::mutate(., group1 = as.character(group1))%>%
  dplyr::mutate(., group2 = as.character(group2))
EC50_stats <- aov(EC50 ~ strain, data = collectEC)%>%
  rstatix::tukey_hsd()%>%
  dplyr::mutate(p.adj.signif = case_when(p.adj > 0.05 ~ "ns",
                                         (p.adj < 0.05 & p.adj > 0.01) ~"*",
                                         (p.adj < 0.01 & p.adj > 0.001) ~"**",
                                         (p.adj < 0.001 & p.adj > 0.0001) ~"***",
                                         (p.adj < 0.0001 & p.adj > 0.00001) ~"****",
                                         (p.adj < 0.00001) ~ "*****"))%>%
  dplyr::filter(grepl("N2",group1) | grepl("N2",group2))

EC50plot <-collectEC%>%
  dplyr::filter(!is.na(strain))%>%
  dplyr::filter(!is.na(EC50))%>%
  dplyr::filter(strain %in% plate_1strains)%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "882","919","920","1081","1082","1075","1076","1327","1328", "1325","1326")))%>%
  ggplot()+
  aes(x=fancy_strain, y = EC50)+
  geom_jitter(width = 0.1,size=0.1)+
  stat_pvalue_manual(EC50_stats, label = "p.adj.signif",xmax="group1", y.position = c(90),remove.bracket = TRUE)+
  geom_boxplot(aes(fill=strain, alpha=0.1), outlier.shape = NA)+
  stat_pvalue_manual(allele_F, label = "strain", xmax = "group2", y.position = c(95), remove.bracket = TRUE)+
  theme_cowplot(8)+
  ylab(glue('EC50 µM FBZ'))+
  scale_x_discrete(labels=c("N2" = "N2", "882" = expression(paste(italic("ean64"))),"919"=expression(paste(italic("ean100"))),"920" = expression(italic("ean101")),"1081"=expression(italic("ean149")),"1082"=expression(italic("ean150")),"1327"=expression(italic("ean203")),"1328"=expression(italic("ean204")),"1325"=expression(italic("ean201")),"1326"=expression(italic("ean202")),"1075"=expression(italic("ean143")),"1076" = expression(italic("ean144"))))+
  scale_fill_manual(name = "fancy_strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V"), values = c("N2" = "orange","882"="grey","919"="red","920"="red","1325"="yellow","1326"="yellow","1327"="purple","1328"="purple","1082" = "blue","1081" = "blue", "1076" = "green","1075" = "green"))+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(size=12,face = "bold",angle = 45,hjust = 1,vjust = 1,margin = unit(c(0,0,0,0),units = "in")),
        axis.text.y = element_text(size=12, face = "bold",margin = unit(c(0 ,0,0,0),units = "in")),
        axis.title.y=element_text(size = 12,face = "bold",margin = unit(c(0,0.075,0,0),units = "in")),
        legend.key.size = unit(0.4,"in"),
        axis.text.y.right = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y.right = element_text(margin = unit(c(0,0,0,0),units = "in"), hjust = 0.1),
        axis.ticks.y.right = element_blank(),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        legend.position = "None")
ggsave(filename = "~/Desktop/ben1_2020_CMD/manuscript/Supplemental_figure_5.jpeg", plot = EC50plot,device = "jpeg",width = 7.5,height = 4,units = "in")
