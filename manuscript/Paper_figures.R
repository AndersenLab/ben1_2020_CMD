#Source the theme 
library(tidyverse)
library(drc)
library(ggpubr)
library(glue)
source("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Scripts/theme.R")
source("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Scripts/scripts.R")
parafull <- c("N2","882","919","920","1075","1076","1081","1082","1325","1326","1327","1328")
#Choose colors for strains 
cols <- c("N2" = "orange","882"="grey","919"="red","920"="red","1325"="yellow","1327"="purple","1082" = "blue","1081" = "blue", "1076" = "green","1075" = "green")
para <- c("N2","882","919","1082","1076","1325","1327")
plate_1strains <- c("N2","882","919","920","1081","1082","1075","1076","1325","1326","1327","1328")
traita <- c("median.EXT")
standerror <- function(x) sd(x)/sqrt(length(x))

##Figure 1
###A
####Load ABZ V3 dose
load("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Processed_data/20190617subtract_pruned.Rdata")
abz_outpruned <- x
rm(x)
abz_outprunedmut <- abz_outpruned%>%
  dplyr::filter(condition != "200")%>%
  dplyr::filter(trait == traita)%>%
  dplyr::mutate(numcon = case_when(condition == "DMSO" ~ 0,
                                   condition == "6_25" ~ 6.25,
                                   condition == "12_5" ~ 12.5,
                                   condition == "25" ~ 25,
                                   condition == "50" ~ 50,
                                   condition == "100" ~ 100))%>%
  dplyr::group_by(strain, condition)%>%
  dplyr::mutate(meancondition = mean(phenotype))%>%
  dplyr::mutate(SND = sd(phenotype))

abzmedian.EXT_subtracted <- abz_outprunedmut%>%
  dplyr::filter(!is.na(strain))%>%
  dplyr::filter(strain %in% para)%>%
  dplyr::filter(condition != "200")%>%
  dplyr::filter(trait == traita)%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "882","919","920","1081","1082","1327","1328", "1325","1326","1075","1076")))%>%
  dplyr::mutate(con_fix = factor(condition, levels = c("DMSO","6_25","12_5","25","50","100")))%>%
  ggplot()+
  geom_errorbar(aes(color = fancy_strain, ymin=meancondition - SND, ymax=meancondition + SND), width=3, size = 0.25)+
  aes(x=numcon, y = phenotype)+
  geom_line(aes(x = numcon, y = meancondition, colour= fancy_strain), size = 0.5)+
  theme_cowplot(8)+
  scale_y_continuous(limits = c(-600,100),sec.axis = dup_axis(name = "Albendazole"))+
  ylab("Normalized animal length")+
  scale_x_continuous(breaks = c(0,6.25,12.5,25,50,100),expand = c(0,0.5))+
  scale_fill_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V"), values = cols)+
  scale_color_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V"), values = cols)+
  base_theme+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, face = "bold",margin = unit(c(0 ,0,0,0),units = "in")),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size = 0.5),
        axis.title.y=element_text(size=12,face = "bold",vjust=-3,margin = unit(c(0,0.075,0,0),units = "in")),
        axis.title.y.right=element_text(size=12,face = "bold",margin = unit(c(0,0,0,0.01),units = "in")),
        legend.key.size = unit(1,"in"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        plot.margin = unit(c(0, 0, 0, 0), "in"),
        legend.position = "None")

###B
####Load FBZ V3
para <- c("N2","882","919","1082","1076","1325","1327")
load("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Processed_data/20190618subtract_pruned.Rdata")
fbz_outpruned <- x

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
  dplyr::mutate(meancondition = mean(phenotype))%>%
  dplyr::mutate(SND = sd(phenotype))

FBZmedian.EXT_subtracted <-fbz_outprunedmut%>%
  dplyr::filter(!is.na(strain))%>%
  dplyr::filter(strain %in% para)%>%
  #dplyr::filter(strain == "N2" | strain == "919" | strain == "1325" | strain == "1327")%>%
  dplyr::filter(condition != "200")%>%
  dplyr::filter(trait == traita)%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "882","919","920","1081","1082","1327","1328", "1325","1326","1075","1076")))%>%
  dplyr::mutate(con_fix = factor(condition, levels = c("DMSO","6_25","12_5","25","50","100")))%>%
  ggplot()+
  geom_errorbar(aes(color = fancy_strain, ymin=meancondition - SND, ymax=meancondition + SND), width=3, size = 0.25)+
  aes(x=numcon, y = phenotype)+
  geom_line(aes(x = numcon, y = meancondition, colour= fancy_strain), size = 0.5)+
  theme_cowplot(8)+
  scale_y_continuous(limits = c(-800,100),sec.axis = dup_axis(name = "Fenbendazole"))+
  ylab("Normalized animal length")+
  xlab("Benzimidazole concentration (µM)")+
  scale_x_continuous(breaks = c(0,6.25,12.5,25,50,100), labels = c("0","6.25","12.5","25","50","100"),expand = c(0,0.5))+
  scale_fill_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V"), values = cols)+
  scale_color_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V"), values = cols)+
  base_theme+
  theme(axis.title.x=element_text(face = "bold", size = 12, margin = unit(c(0,0,0,0), units = "in"),vjust=1),
        axis.text.x = element_text(face = "bold", margin = unit(c(0.05,0,0,0), units = "in"), angle = 45),
        axis.text.y = element_text(size=12, face = "bold",margin = unit(c(0 ,0,0,0),units = "in")),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.ticks.x = element_line(size = 0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.title.y=element_text(size=12,face = "bold",vjust=-3,margin = unit(c(0,0.075,0,0),units = "in")),
        axis.title.y.right=element_text(size=12,face = "bold", margin = unit(c(0,0,0,0.01),units = "in")),
        legend.key.size = unit(0.4,"in"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        plot.margin = unit(c(0, 0, 0, 0), "in"), 
        legend.position = "None")

##Plot and save figure
doseresponses <- plot_grid(abzmedian.EXT_subtracted,FBZmedian.EXT_subtracted,nrow = 2,labels = c("A","B"),vjust = 1,hjust = 0.075,axis = "l",align="v",label_size = 12,label_fontfamily = "arial")

ggsave("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/manuscript/Figure1.png",plot = doseresponses,device = "png",width = 3.75,height = 4,units = "in")

##Figure 2
load("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Processed_data/20190618subtract_pruned.Rdata")
fbzfordoseresponse <-x%>%
  dplyr::filter(!is.na(strain))%>%
  #dplyr::filter(strain %in% para)%>%
  dplyr::filter(condition != "200")%>%
  dplyr::filter(trait == traita)%>%
  dplyr::mutate(cond = case_when(condition == "DMSO" ~ 0,
                                 condition == "6_25" ~ 6.25,
                                 condition == "12_5" ~ 12.5,
                                 condition == "25" ~ 25,
                                 condition == "50" ~ 50,
                                 condition == "100" ~ 100))%>%
  dplyr::mutate(strain_comb = case_when(strain == "N2" ~ "N2",
                                        strain == "882" ~ "del",
                                        strain == "919" ~ "F200Y",
                                        strain == "1076"  ~ "F167Y",
                                        strain == "1082" ~ "E198A",
                                        strain == "1325" ~ "E198L",
                                        strain == "1327" ~ "E198V",
                                        strain == "920" ~ "920",
                                        strain == "1075" ~ "1075",
                                        strain == "1081" ~ "1081",
                                        strain == "1326" ~ "1326",
                                        strain == "1328" ~ "1328"))%>%
  dplyr::mutate(medianEXT_norm = ((phenotype - min(phenotype)) / (max(phenotype) - min(phenotype))))
fbzfull <- NULL
for(i in unique(fbz_outprunedmut$strain)){
  try <- fbzfordoseresponse %>%
    dplyr::filter(strain == i)%>%
    dplyr::mutate( q90_norm = (phenotype - min(phenotype)) / (max(phenotype) - min(phenotype)))
  fbzfull <- rbind(fbzfull,try)
  rm(try)
}
rm(x)
fbzec50 <- drc::drm(q90_norm ~ cond,strain_comb, data = fbzfull, fct = LL.3())


summary(fbzec50)
ec50 <- ED(fbzec50,50)
ec50fix <- data.table::setDT(data.frame(ec50), keep.rownames = TRUE)%>%
  dplyr::mutate(strain = factor(rn, levels = c("e:N2:50","e:del:50","e:F200Y:50","e:E198A:50","e:F167Y:50","e:E198V:50","e:E198L:50")))


FBZec50 <- ec50fix %>%
  dplyr::filter(!is.na(strain))%>%
  ggplot()+
  aes(x = strain, y = Estimate)+
  geom_col(aes(fill = strain))+
  #stat_compare_means(paired = TRUE,method = "t.test")+
  geom_errorbar(aes(ymin=Estimate - Std..Error, ymax = Estimate + Std..Error), width = 0.2)+
  scale_x_discrete(labels=c("e:N2:50" = "N2", "e:del:50" = "del","e:F200Y:50" = "F200Y","e:E198A:50"="E198A","e:F167Y:50"="F167Y","e:E198V:50"="E198V","e:E198L:50"="E198L"))+
  scale_fill_manual( values = c("e:N2:50"="orange","e:del:50"="grey","e:F200Y:50"="red","e:E198A:50"="blue","e:F167Y:50"="green","e:E198V:50"="purple","e:E198L:50"="yellow"))+
  #scale_y_continuous(limits = c(20,45), oob = rescale_none)+
  xlab("Strain")+
  ylab("Estimated EC50 (µM FBZ)")+
  base_theme+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(face = "bold", angle = 45,hjust = 1,vjust = 1,margin = unit(c(0,0,0,0),units = "in")),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.text.y=element_text(size=12, face = "bold",margin = unit(c(0 ,0,0,0),units = "in")),
        axis.title.y=element_text(size = 12, face = "bold",margin = unit(c(0,0.075,0,0),units = "in")),
        legend.key.size = unit(0.4,"in"),
        legend.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.text = element_text(face = "bold"),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "in"),
        legend.position = "None")
ggsave("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/manuscript/Figure2.png",FBZec50, device = "png",width = 3.75,height = 2,units = "in")  

   ##Figure 3
###A
####Load DMSO from abz data for high throughput 
load("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Processed_data/20200130_noassayregressed.Rdata")
DMSOabz <- biopruned %>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(strain %in% para)
rm(biopruned)

dmsoabz_stats <- aov(median.EXT ~ strain, data = DMSOabz)%>%
  rstatix::tukey_hsd()%>%
  dplyr::filter(grepl("N2",group1) | grepl("N2",group2))

DMSOabz_medianEXT <-DMSOabz%>%
  dplyr::filter(!is.na(strain))%>%
  dplyr::filter(strain %in% para)%>%
  dplyr::mutate(medianEXT = median.EXT/10)%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "882","919","920","1081","1082","1075","1076","1327","1328", "1325","1326")))%>%
  ggplot()+
  aes(x=fancy_strain, y = median.EXT)+
  geom_jitter(width = 0.1,size=0.1)+
  #geom_violin(aes(fill=fancy_strain,alpha=0.1))+
  geom_boxplot(aes(fill=strain,alpha=0.1), outlier.shape = NA)+
  theme_cowplot(8)+
  stat_pvalue_manual(dmsoabz_stats, label = "p.adj.signif",xmax="group1", y.position = c(1525),remove.bracket = TRUE)+
  ylab(glue('Animal length'))+
  scale_y_continuous(limits = c(900,1530),sec.axis = dup_axis(name = "DMSO"))+
  scale_x_discrete(labels=c("N2" = "N2", "882" = "del","919" = "F200Y","1082"="E198A","1327"="E198V","1325"="E198L","1076"="F167Y"))+
  scale_fill_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V"), values = cols)+
  scale_color_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V"), values = cols)+
  base_theme+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=12, face = "bold",margin = unit(c(0 ,0,0,0),units = "in")),
        axis.title.y=element_text(size = 12,face = "bold",vjust = -2,margin = unit(c(0,0.075,0,0),units = "in")),
        legend.key.size = unit(0.4,"in"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(margin = unit(c(0,0,0,0),units = "in")),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        plot.margin = unit(c(0.1, 0.1, 0, 0), "in"),
        legend.position = "None")
###B
####Load DMSO from fbz data for high throughput
load("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Processed_data/20190606pruned_constrain.Rdata")
plates1 <-c("1","2","3","4","5","6","7","8","9","10")
parafbz <- c("N2","882","920","1081","1075","1325","1327")
DMSOfbz <- out_pruned %>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::filter(trait == traita)%>%
  dplyr::filter(strain %in% parafbz)%>%
  dplyr::filter(plate %in% plates1)

dmsofbz_stats <- aov(phenotype ~ strain, data = DMSOfbz)%>%
  rstatix::tukey_hsd()%>%
  dplyr::filter(grepl("N2",group1) | grepl("N2",group2))
rm(out_pruned)

DMSOfbz_medianEXT <-DMSOfbz%>%
  dplyr::filter(!is.na(strain))%>%
  dplyr::filter(trait == traita)%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "882","919","920","1081","1082","1075","1076","1327","1328", "1325","1326")))%>%
  ggplot()+
  aes(x=fancy_strain, y = phenotype)+
  geom_jitter(width = 0.1,size=1)+
  geom_boxplot(aes(fill=strain,alpha=0.1),outlier.shape = NA)+
  theme_cowplot(8)+
  ylim(1005,1300)+
  scale_y_continuous(sec.axis = dup_axis(name = "DMSO"))+
  ylab(glue('{traita} [DMSO]'))+
  stat_pvalue_manual(dmsofbz_stats, label = "p.adj.signif",xmax="group1", y.position = c(1275),remove.bracket = TRUE)+
  scale_x_discrete(labels=c("N2" = "N2", "882" = "del","920" = "F200Y","1081"="E198A","1327"="E198V","1325"="E198L","1075"="F167Y"))+
  scale_fill_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V"), values = cols)+
  scale_color_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V"), values = cols)+
  base_theme+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(face = "bold"),
        axis.ticks.x = element_blank(),
        axis.text.y=element_text(size=12, face = "bold",margin = unit(c(0 ,0,0,0),units = "in")),
        axis.title.y=element_text(size = 12,face = "bold",vjust = -2,margin = unit(c(0,0.075,0,0),units = "in")),
        legend.key.size = unit(0.4,"in"),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(margin = unit(c(0,0,0,0),units = "in")),
        axis.ticks.y.right = element_blank(),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        plot.margin = unit(c(0.1, 0.1, 0, 0), "in"),
        legend.position = "None")
###C
####Load ABZ data for high throughput
load("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Processed_data/20200130.Rdata")
preregressedabz <- out_pruned %>%
  dplyr::filter(strain %in% para)%>%
  dplyr::filter(trait == traita)
regressedabz <- preregressedabz %>%
  dplyr::mutate(assay = "A")%>%
  easysorter::regress(.,assay = FALSE)
rm(out_pruned)

abz_stats <- aov(phenotype ~ strain, data = regressedabz)%>%
  rstatix::tukey_hsd()%>%
  dplyr::filter(grepl("N2",group1) | grepl("N2",group2))


regressedabz_medianEXT <-regressedabz%>%
  dplyr::filter(!is.na(strain))%>%
  dplyr::filter(trait == traita)%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "882","919","920","1081","1082","1075","1076","1327","1328", "1325","1326")))%>%
  ggplot()+
  aes(x=fancy_strain, y = phenotype)+
  #geom_violin(aes(fill=fancy_strain,alpha=0.1))+
  geom_jitter(width = 0.1,size=0.1)+
  geom_boxplot(aes(fill=strain, alpha = 0.1), outlier.shape = NA)+
  theme_cowplot(8)+
  scale_y_continuous(sec.axis = dup_axis(name = "Albendazole"))+
  stat_pvalue_manual(abz_stats, label = "p.adj.signif",xmax="group1", y.position = c(275),remove.bracket = TRUE)+
  ylab(glue("Normalized animal length"))+
  scale_x_discrete(labels=c("N2" = "N2", "882" = "del","919" = "F200Y","1082"="E198A","1327"="E198V","1325"="E198L","1076"="F167Y"))+
  scale_fill_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V"), values = cols)+
  scale_color_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V"), values = cols)+
  base_theme+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=12, face = "bold",margin = unit(c(0 ,0,0,0),units = "in")),
        axis.title.y=element_text(size = 12,face = "bold",vjust = -2,margin = unit(c(0,0.075,0,0),units = "in")),
        legend.key.size = unit(0.4,"in"),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(margin = unit(c(0,0,0,0),units = "in")),
        axis.ticks.y.right = element_blank(),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        plot.margin = unit(c(0.1, 0.1, 0, 0), "in"),
        legend.position = "None")
###D
####Load FBZ data for high throughput
load("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Processed_data/20190606pruned_constrain.Rdata")
preregressedfbz <- out_pruned %>%
  dplyr::filter(strain %in% parafbz)%>%
  dplyr::filter(plate %in% plates1)

regressedfbz <- easysorter::regress(preregressedfbz)

fbz <- regressedfbz%>%
  dplyr::filter(trait == traita)%>%
  dplyr::filter(strain %in% parafbz)

fbz_stats <- aov(phenotype ~ strain, data = fbz)%>%
  rstatix::tukey_hsd()%>%
  dplyr::filter(grepl("N2",group1) | grepl("N2",group2))


regressedfbz_medianEXT <-regressedfbz%>%
  dplyr::filter(!is.na(strain))%>%
  dplyr::filter(trait == traita)%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "882","919","920","1081","1082","1075","1076","1327","1328", "1325","1326")))%>%
  ggplot()+
  aes(x=fancy_strain, y = phenotype)+
  geom_jitter(width = 0.1,size=0.1)+
  #geom_violin(aes(fill=fancy_strain,alpha=0.1))+
  geom_boxplot(aes(fill=strain, alpha=0.1), outlier.shape = NA)+
  theme_cowplot(8)+
  ylab(glue('Normalized animal length'))+
  scale_y_continuous(limits = c(-275,200),sec.axis = dup_axis(name = "Fenbendazole"))+
  stat_pvalue_manual(fbz_stats, label = "p.adj.signif",xmax="group1", y.position = c(190),remove.bracket = TRUE)+
  scale_x_discrete(labels=c("N2" = "N2", "882" = "del","920" = "F200Y","1081"="E198A","1327"="E198V","1325"="E198L","1075"="F167Y"))+
  scale_fill_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V"), values = cols)+
  scale_color_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V"), values = cols)+
  base_theme+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(face = "bold",angle = 45,hjust = 1,vjust = 1,margin = unit(c(0,0,0,0),units = "in")),
        axis.text.y = element_text(size=12, face = "bold",margin = unit(c(0 ,0,0,0),units = "in")),
        axis.title.y=element_text(size = 12,face = "bold",vjust = -2,margin = unit(c(0,0.075,0,0),units = "in")),
        legend.key.size = unit(0.4,"in"),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(margin = unit(c(0,0,0,0),units = "in")),
        axis.ticks.y.right = element_blank(),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0), "in"),
        legend.position = "None")

##Plot and save figure

highthroughput <- plot_grid(DMSOabz_medianEXT,regressedabz_medianEXT,regressedfbz_medianEXT,nrow = 3,ncol = 1,labels=c("A","B","C"),vjust = 1,hjust = 0.075,axis = "l",align="v", label_size = 12,label_fontfamily = "arial",rel_widths = c(1,1,1),rel_heights = c(1,1,1))

ggsave("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/manuscript/Figure3.png",plot = highthroughput,device = "png",width = 3.75,height = 6,units = "in")


##Figure 4
###A
####Load DMSO comp data
comp_data <- readr::read_tsv("~/Dropbox/AndersenLab/LabFolders/Clay/2019_ABZcomp/untitled folder/compdata_full.tsv")
comp_data$Allele_Freq <- as.numeric(comp_data$Allele_Freq)
comp_data$Mean <- as.numeric(comp_data$Mean)
cols <- c("N2" = "orange", "del" = "grey", "F200Y" = "Red", "E198V"= "Purple","E198L"="yellow","E198A"="blue","F167Y"="green")
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
  geom_errorbar(aes(ymin=Mean - SND, ymax=Mean + SND), width=0.5,size=0.5)+
  geom_line(aes(linetype = Condition_f),size=1)+
  geom_point()+
  scale_color_manual(labels = c("N2"="N2","del" = "del","F200Y"="F200Y","E198V"="E198V","E198L"="E198L","E198A" = "E198A","F167Y"="F167Y"), values = cols)+
  theme_cowplot(8)+
  scale_y_continuous(breaks = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),limits = c(0.3,1.0))+
  base_theme+
  ylab("Relative allele frequency")+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(face = "bold",margin = unit(c(0,0,0,0),units = "in")),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.text.y=element_text(size=12, face = "bold",margin = unit(c(0 ,0,0,0),units = "in")),
        axis.title.y=element_text(size = 12, face = "bold",margin = unit(c(0,0.075,0,0),units = "in")),
        legend.key.size = unit(0.4,"in"),
        plot.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.text = element_text(face = "bold"),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "in"),
        legend.position = "None")
###B
####Load DMSO box plots
fitness <- readr::read_tsv("~/Dropbox/AndersenLab/LabFolders/Clay/2019_ABZcomp/untitled folder/fitness_fixed_full.txt")
Dfitness <- dplyr::filter(fitness,Condition == "D")
dmso_stats <- aov(Fitness ~ Strain, data = Dfitness)%>%
  rstatix::tukey_hsd()%>%
  dplyr::filter(grepl("N2",group1) | grepl("N2",group2))
DFIT <- fitness %>%
  dplyr::filter(Condition == "D")%>%
  dplyr::mutate(fancy_strain=factor(Strain, levels = c("N2", "del","F200Y","E198A","F167Y","E198V","E198L")))%>%
  ggplot()+
  aes(x=fancy_strain, y=Fitness)+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_jitter(width = 0.1)+
  geom_boxplot(aes(fill = Strain, alpha=0.1),outlier.shape = NA)+
  theme_cowplot(8)+
  stat_pvalue_manual(dmso_stats, label = "p.adj.signif",xmax="group1", y.position = c(0.15),remove.bracket = TRUE)+
  base_theme+
  ylab("Fitness")+
  scale_x_discrete(labels = c("N2"="N2","del" = "del","F200Y"="F200Y","E198A" ="E198A","F167Y"="F167Y","E198L"="E198L","E198V"="E198V"))+
  scale_color_manual(labels = c("N2"="N2","del" = "del","F200Y"="F200Y","E198V"="E198V"), values = cols)+
  scale_fill_manual(labels = c("N2"="N2","del" = "del","F200Y"="F200Y","E198V"="E198V"), values = cols)+
  theme(legend.position="none")+
  ylab("Fitness")+
  xlab("Strain")+
  scale_y_continuous(sec.axis = dup_axis(name = "DMSO"))+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(face = "bold", angle = 45,hjust = 1,vjust = 1,margin = unit(c(0,0,0,0),units = "in")),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.text.y=element_text(size=12, face = "bold",margin = unit(c(0 ,0,0,0),units = "in")),
        axis.title.y=element_text(size = 12, face = "bold",margin = unit(c(0,0.075,0,0),units = "in")),
        legend.key.size = unit(0.4,"in"),
        legend.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.text = element_text(face = "bold"),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "in"),
        legend.position = "None")

####Load ABZ comp data
ABZ_comp <- x %>%
  dplyr::filter(Condition_f == "Albendazole")%>%
  ggplot()+
  aes(x = Generation, y = Mean , fill = Condition_f, color = Strain) +
  geom_errorbar(aes(ymin=Mean - SND, ymax=Mean + SND), width=0.5,size=0.5)+
  geom_line(aes(linetype = Condition_f), size=1)+
  geom_point()+
  scale_color_manual(labels = c("N2"="N2","del" = "del","F200Y"="F200Y","E198V"="E198V","E198L"="E198L","E198A" = "E198A","F167Y"="F167Y"), values = cols)+
  theme_cowplot(8)+
  scale_y_continuous(breaks = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),limits = c(0.3,1.0))+
  base_theme+
  ylab("Relative allele frequency")+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(face = "bold",margin = unit(c(0,0,0,0),units = "in")),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.text.y=element_text(size=12, face = "bold",margin = unit(c(0 ,0,0,0),units = "in")),
        axis.title.y=element_text(size = 12, face = "bold",margin = unit(c(0,0.075,0,0),units = "in")),
        legend.key.size = unit(0.4,"in"),
        plot.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.text = element_text(face = "bold"),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "in"),
        legend.position = "None")
###D
####Load ABZ box plots
Afitness <- dplyr::filter(fitness,Condition == "A")
A_stats <- aov(Fitness ~ Strain, data = Afitness)%>%
  rstatix::tukey_hsd()%>%
  dplyr::filter(grepl("N2",group1) | grepl("N2",group2))
A_F200YE198Vstats <- aov(Fitness ~ Strain, data = Afitness)%>%
  rstatix::tukey_hsd()%>%
  dplyr::filter(grepl("E198V",group1) , grepl("F200Y",group2))
AFIT <- fitness %>%
  dplyr::filter(Condition == "A")%>%
  dplyr::mutate(fancy_strain=factor(Strain, levels = c("N2", "del","F200Y","E198A","F167Y","E198V","E198L")))%>%
  ggplot()+
  aes(x=fancy_strain, y=Fitness)+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_jitter(width = 0.1)+
  geom_boxplot(aes(fill = Strain, alpha=0.1),outlier.shape = NA)+
  theme_cowplot(8)+
  stat_pvalue_manual(A_stats, label = "p.adj.signif",xmax="group1", y.position = c(0.75),remove.bracket = TRUE)+
  stat_pvalue_manual(A_F200YE198Vstats, label = "p.adj.signif",y.position = c(0.70), remove.bracket = FALSE)+
  base_theme+
  scale_x_discrete(labels = c("N2"="N2","del" = "del","F200Y"="F200Y","E198A" ="E198A","F167Y"="F167Y","E198L"="E198L","E198V"="E198V"))+
  scale_color_manual(labels = c("N2"="N2","del" = "del","F200Y"="F200Y","E198V"="E198V"), values = cols)+
  scale_fill_manual(labels = c("N2"="N2","del" = "del","F200Y"="F200Y","E198V"="E198V"), values = cols)+
  theme(legend.position="none")+
  scale_y_continuous(sec.axis = dup_axis(name = "Albendazole"))+
  ylab("Fitness")+
  xlab("Strain")+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(face = "bold", angle = 45,hjust = 1,vjust = 1,margin = unit(c(0,0,0,0),units = "in")),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.text.y=element_text(size=12, face = "bold",margin = unit(c(0 ,0,0,0),units = "in")),
        axis.title.y=element_text(size = 12, face = "bold",margin = unit(c(0,0.075,0,0),units = "in")),
        legend.key.size = unit(0.4,"in"),
        legend.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.text = element_text(face = "bold"),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "in"),
        legend.position = "None")
##Plot and save plot
competitionplot <- plot_grid(DMSO_comp,DFIT,ABZ_comp,AFIT,nrow = 2,ncol = 2,labels = c("A","B","C","D"),align = "vh",axis = "lrbt",label_size = 12,label_fontfamily = "arial",rel_widths = c(1,1,1,1),rel_heights = c(1,1,1,1))
ggsave("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/manuscript/Figure4.png",plot = competitionplot,device = "png",width = 7.5,height = 6,units = "in")
##Figure 5
###A


