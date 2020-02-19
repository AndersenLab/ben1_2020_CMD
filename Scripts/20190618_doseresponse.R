#Theme Options
# colors
axis_color <- "#000F08"
highlight_color <- "#D7263D"
background_color <- "white"

# font
nµMber_font <- "Helvetica"
axes_text_size <- 20
axes_title_font <- "Helvetica"
axes_title_size <- 20
title_size <- 20

base_theme <- theme(
  line = element_line(colour = axis_color, size = 0.5, linetype = 1), 
  rect = element_rect(fill = background_color, colour = axis_color, size = 0.5, linetype = 1), 
  text = element_text(family = axes_title_font, size = axes_text_size), 
  
  axis.text = element_text(family = nµMber_font,size = rel(0.8), colour = "grey30", margin = unit(0.1, "cm")),
  axis.text.x = element_text(vjust = 1), 
  axis.text.y = element_text(hjust = 1), 
  axis.ticks = element_line(colour = "gray90"), 
  axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), 
  axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), angle = 90), 
  axis.ticks.length = unit(0.15, "cm"),
  
  strip.text = element_text(size = rel(0.8)), 
  strip.background = element_rect(fill = background_color, colour = NA, size = 0.5, linetype = 1),
  
  plot.background = element_rect(fill = background_color, color = NA),
  
  legend.background = element_rect(fill=background_color,colour = background_color), 
  legend.spacing = unit(0.2, "cm"), 
  legend.key = element_rect(fill = background_color, colour = NA), 
  legend.key.size = unit(1, "lines"), 
  legend.key.height = NULL, 
  legend.key.width = NULL, 
  legend.text = element_text(size = rel(0.6)), 
  legend.text.align = NULL, 
  legend.title = element_text(size = rel(0.6), hjust = 0), 
  legend.title.align = NULL, 
  legend.position = "right", 
  legend.direction = NULL, 
  legend.justification = "center", 
  legend.box = NULL, 
  
  panel.background = element_rect(fill = background_color, colour = NA),  
  panel.grid.major = element_line(colour = "gray90"), 
  panel.grid.minor = element_blank(), 
  panel.spacing = unit(1, "lines"), 
  panel.margin.x = NULL, 
  panel.margin.y = NULL)

load("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Processed_data/20190618subtract_constrain.Rdata")
#Pick Colors for strains
cols <- c("N2" = "orange","882"="grey","919"="red", "920" = "red", "1081" = "blue", "1075" = "green", "1325" = "brown","1326"="brown", "1327" = "yellow","N2" = "orange", "919" = "red", "1082" = "blue", "1076" = "green", "1326" = "brown", "1328" = "yellow","1139"="black","1140"="black","1137"="pink","1138"="pink","1314"="tan","1317"="tan","1097"="maroon","1098"="maroon")
#Parasite strains
p <- c("N2", "882","919","920", "1325","1326","1327","1328","1081","1082","1076","1075")

#Plot normalized difference between N2 and the deletion strain
r01_strains <- c("N2","882")
first <- subtract_dr_to_pc%>%
  dplyr::filter(strain %in% r01_strains, !is.na(q90.EXT))
test <- first%>%
  dplyr::select(strain,condition, q90.EXT)%>%
  dplyr::filter(!is.na(q90.EXT))%>%
  dplyr::mutate(mednorm = ((first$q90.EXT - min(first$q90.EXT)) / (max(first$q90.EXT) - min(first$q90.EXT))))%>%
  dplyr::filter(mednorm < 0.8)
  
test2 <- test%>%
  dplyr::select(strain,condition, q90.EXT, mednorm)%>%
  dplyr::mutate(mednorm2 = ((test$mednorm - min(test$mednorm)) / (max(test$mednorm) - min(test$mednorm))))


test2%>%
  filter(condition!="200",
         #median.EXT < 500,
         #median.EXT >-1000,
         strain %in% r01_strains,
         #strain %in% plate_set_2,
         #!(plate %in% plates1),
         !is.na(strain))%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "882","919","920","1081","1082","1327","1328", "1325","1326","1075","1076","1137","1138","1139","1140","1314","1317","1097","1098")))%>%
  dplyr::mutate(con_fix = factor(condition, levels = c("DMSO","6_25","12_5","25","50","100","200")))%>%
  ggplot()+
  aes(x=con_fix, y = mednorm2)+
  geom_smooth(alpha=0.1, se=TRUE, aes(group=fancy_strain, colour=fancy_strain), method = 'auto')+
  theme_bw(24)+
  #geom_jitter()+
  ylab("median EXT")+
  scale_x_discrete(labels=c("DMSO" = "DMSO", "6_25" = "6.25µM","12_5" = "12.5µM","25"="25µM","50"="50µM","100"="100µM"))+
  scale_fill_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y","1138"="A185P","1140"="A185P","1138"="Q131L","1098"="S145F", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V","1137"="Q131L","1139"="A185P","1314"="M257I","1317"="M257I","1097"="S145F"), values = cols)+
  scale_color_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y","1138"="A185P","1140"="A185P","1138"="Q131L","1098"="S145F", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V","1137"="Q131L","1139"="A185P","1314"="M257I","1317"="M257I","1097"="S145F"), values = cols)+
  base_theme+
  theme(axis.title.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size=28,face = "bold"),
        axis.text.y=element_text(size=28,face = "bold"),
        axis.title.y=element_text(size=30,face = "bold"),
        legend.key.size = unit(0.4,"in"),
        legend.title = element_text(size=28,face = "bold"),
        legend.text = element_text(size=28,face = "bold"),
        legend.position = "None")

#Plot DMSO conditions Median.EXT
subtract_dr_to_pc%>%
  filter(condition!="200",
         median.EXT < 500,
         median.EXT >-1000,
         strain %in% p,
         #strain %in% plate_set_2,
         #!(plate %in% plates1),
         !is.na(strain))%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "882","919","920","1081","1082","1327","1328", "1325","1326","1075","1076","1137","1138","1139","1140","1314","1317","1097","1098")))%>%
  dplyr::mutate(con_fix = factor(condition, levels = c("DMSO","6_25","12_5","25","50","100","200")))%>%
  ggplot()+
  aes(x=con_fix, y = median.EXT)+
  geom_smooth(alpha=0.1, se=TRUE, aes(group=fancy_strain, colour=fancy_strain), method = 'auto')+
  theme_bw(24)+
  geom_jitter()+
  ylab("median EXT")+
  scale_x_discrete(labels=c("DMSO" = "DMSO", "6_25" = "6.25µM","12_5" = "12.5µM","25"="25µM","50"="50µM","100"="100µM"))+
  scale_fill_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y","1138"="A185P","1140"="A185P","1138"="Q131L","1098"="S145F", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V","1137"="Q131L","1139"="A185P","1314"="M257I","1317"="M257I","1097"="S145F"), values = cols)+
  scale_color_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y","1138"="A185P","1140"="A185P","1138"="Q131L","1098"="S145F", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V","1137"="Q131L","1139"="A185P","1314"="M257I","1317"="M257I","1097"="S145F"), values = cols)+
  base_theme+
  theme(axis.title.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size=28,face = "bold"),
        axis.text.y=element_text(size=28,face = "bold"),
        axis.title.y=element_text(size=30,face = "bold"),
        legend.key.size = unit(0.4,"in"),
        legend.title = element_text(size=28,face = "bold"),
        legend.text = element_text(size=28,face = "bold"),
        legend.position = "None")
ggsave("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/plots/20190618_FBZdose/subtract_dose_medianEXT.png", device = "png", plot = last_plot(), width = 14, height = 6, units = "in")

#Plot DMSO conditions Median.TOF
subtract_dr_to_pc%>%
  filter(condition!="200",
         median.EXT < 500,
         median.EXT >-1000,
         strain %in% p,
         #strain %in% plate_set_2,
         #!(plate %in% plates1),
         !is.na(strain))%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "882","919","920","1081","1082","1327","1328", "1325","1326","1075","1076","1137","1138","1139","1140","1314","1317","1097","1098")))%>%
  dplyr::mutate(con_fix = factor(condition, levels = c("DMSO","6_25","12_5","25","50","100","200")))%>%
  ggplot()+
  aes(x=con_fix, y = median.TOF)+
  geom_smooth(alpha=0.1, se=TRUE, aes(group=fancy_strain, colour=fancy_strain), method = 'auto')+
  theme_bw(24)+
  #geom_jitter()+
  ylab("median TOF")+
  scale_x_discrete(labels=c("DMSO" = "DMSO", "6_25" = "6.25µM","12_5" = "12.5µM","25"="25µM","50"="50µM","100"="100µM"))+
  scale_fill_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y","1138"="A185P","1140"="A185P","1138"="Q131L","1098"="S145F", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V","1137"="Q131L","1139"="A185P","1314"="M257I","1317"="M257I","1097"="S145F"), values = cols)+
  scale_color_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y","1138"="A185P","1140"="A185P","1138"="Q131L","1098"="S145F", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V","1137"="Q131L","1139"="A185P","1314"="M257I","1317"="M257I","1097"="S145F"), values = cols)+
  base_theme+
  theme(axis.title.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size=28,face = "bold"),
        axis.text.y=element_text(size=28,face = "bold"),
        axis.title.y=element_text(size=30,face = "bold"),
        legend.key.size = unit(0.4,"in"),
        legend.title = element_text(size=28,face = "bold"),
        legend.text = element_text(size=28,face = "bold"),
        legend.position = "None")
ggsave("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/plots/20190618_FBZdose/subtract_dose_medianTOF.png", device = "png", plot = last_plot(), width = 14, height = 6, units = "in")


#Plot DMSO conditions mean.EXT
subtract_dr_to_pc%>%
  filter(condition!="200",
         median.EXT < 500,
         median.EXT >-1000,
         strain %in% p,
         #strain %in% plate_set_2,
         #!(plate %in% plates1),
         !is.na(strain))%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "882","919","920","1081","1082","1327","1328", "1325","1326","1075","1076","1137","1138","1139","1140","1314","1317","1097","1098")))%>%
  dplyr::mutate(con_fix = factor(condition, levels = c("DMSO","6_25","12_5","25","50","100","200")))%>%
  ggplot()+
  aes(x=con_fix, y = mean.EXT)+
  geom_smooth(alpha=0.1, se=TRUE, aes(group=fancy_strain, colour=fancy_strain), method = 'auto')+
  theme_bw(24)+
  #geom_jitter()+
  ylab("mean EXT")+
  scale_x_discrete(labels=c("DMSO" = "DMSO", "6_25" = "6.25µM","12_5" = "12.5µM","25"="25µM","50"="50µM","100"="100µM"))+
  scale_fill_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y","1138"="A185P","1140"="A185P","1138"="Q131L","1098"="S145F", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V","1137"="Q131L","1139"="A185P","1314"="M257I","1317"="M257I","1097"="S145F"), values = cols)+
  scale_color_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y","1138"="A185P","1140"="A185P","1138"="Q131L","1098"="S145F", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V","1137"="Q131L","1139"="A185P","1314"="M257I","1317"="M257I","1097"="S145F"), values = cols)+
  base_theme+
  theme(axis.title.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size=28,face = "bold"),
        axis.text.y=element_text(size=28,face = "bold"),
        axis.title.y=element_text(size=30,face = "bold"),
        legend.key.size = unit(0.4,"in"),
        legend.title = element_text(size=28,face = "bold"),
        legend.text = element_text(size=28,face = "bold"),
        legend.position = "None")
ggsave("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/plots/20190618_FBZdose/subtract_dose_meanEXT.png", device = "png", plot = last_plot(), width = 14, height = 6, units = "in")

#Plot DMSO conditions mean.TOF
subtract_dr_to_pc%>%
  filter(condition!="200",
         median.EXT < 500,
         median.EXT >-1000,
         strain %in% p,
         #strain %in% plate_set_2,
         #!(plate %in% plates1),
         !is.na(strain))%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "882","919","920","1081","1082","1327","1328", "1325","1326","1075","1076","1137","1138","1139","1140","1314","1317","1097","1098")))%>%
  dplyr::mutate(con_fix = factor(condition, levels = c("DMSO","6_25","12_5","25","50","100","200")))%>%
  ggplot()+
  aes(x=con_fix, y = mean.TOF)+
  geom_smooth(alpha=0.1, se=TRUE, aes(group=fancy_strain, colour=fancy_strain), method = 'auto')+
  theme_bw(24)+
  #geom_jitter()+
  ylab("mean TOF")+
  scale_x_discrete(labels=c("DMSO" = "DMSO", "6_25" = "6.25µM","12_5" = "12.5µM","25"="25µM","50"="50µM","100"="100µM"))+
  scale_fill_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y","1138"="A185P","1140"="A185P","1138"="Q131L","1098"="S145F", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V","1137"="Q131L","1139"="A185P","1314"="M257I","1317"="M257I","1097"="S145F"), values = cols)+
  scale_color_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y","1138"="A185P","1140"="A185P","1138"="Q131L","1098"="S145F", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V","1137"="Q131L","1139"="A185P","1314"="M257I","1317"="M257I","1097"="S145F"), values = cols)+
  base_theme+
  theme(axis.title.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size=28,face = "bold"),
        axis.text.y=element_text(size=28,face = "bold"),
        axis.title.y=element_text(size=30,face = "bold"),
        legend.key.size = unit(0.4,"in"),
        legend.title = element_text(size=28,face = "bold"),
        legend.text = element_text(size=28,face = "bold"),
        legend.position = "None")
ggsave("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/plots/20190618_FBZdose/subtract_dose_meanTOF.png", device = "png", plot = last_plot(), width = 14, height = 6, units = "in")



#Plot DMSO conditions q90.tof 
subtract_dr_to_pc%>%
  filter(condition!="200",
         median.EXT < 500,
         median.EXT >-1000,
         strain %in% p,
         #strain %in% plate_set_2,
         #!(plate %in% plates1),
         !is.na(strain))%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "882","919","920","1081","1082","1327","1328", "1325","1326","1075","1076","1137","1138","1139","1140","1314","1317","1097","1098")))%>%
  dplyr::mutate(con_fix = factor(condition, levels = c("DMSO","6_25","12_5","25","50","100","200")))%>%
  ggplot()+
  aes(x=con_fix, y = q90.TOF)+
  geom_smooth(alpha=0.1, se=TRUE, aes(group=fancy_strain, colour=fancy_strain), method = 'auto')+
  theme_bw(24)+
  #geom_jitter()+
  ylab("q90 TOF")+
  scale_x_discrete(labels=c("DMSO" = "DMSO", "6_25" = "6.25µM","12_5" = "12.5µM","25"="25µM","50"="50µM","100"="100µM"))+
  scale_fill_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y","1138"="A185P","1140"="A185P","1138"="Q131L","1098"="S145F", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V","1137"="Q131L","1139"="A185P","1314"="M257I","1317"="M257I","1097"="S145F"), values = cols)+
  scale_color_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y","1138"="A185P","1140"="A185P","1138"="Q131L","1098"="S145F", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V","1137"="Q131L","1139"="A185P","1314"="M257I","1317"="M257I","1097"="S145F"), values = cols)+
  base_theme+
  theme(axis.title.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size=28,face = "bold"),
        axis.text.y=element_text(size=28,face = "bold"),
        axis.title.y=element_text(size=30,face = "bold"),
        legend.key.size = unit(0.4,"in"),
        legend.title = element_text(size=28,face = "bold"),
        legend.text = element_text(size=28,face = "bold"),
        legend.position = "None")
ggsave("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/plots/20190618_FBZdose/subtract_dose_q90tof.png", device = "png", plot = last_plot(), width = 14, height = 6, units = "in")



#Plot DMSO conditions q90.EXT 
subtract_dr_to_pc%>%
  filter(condition!="200",
         median.EXT < 500,
         median.EXT >-1000,
         strain %in% p,
         #strain %in% plate_set_2,
         #!(plate %in% plates1),
         !is.na(strain))%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "882","919","920","1081","1082","1327","1328", "1325","1326","1075","1076","1137","1138","1139","1140","1314","1317","1097","1098")))%>%
  dplyr::mutate(con_fix = factor(condition, levels = c("DMSO","6_25","12_5","25","50","100","200")))%>%
  ggplot()+
  aes(x=con_fix, y = q90.EXT)+
  geom_smooth(alpha=0.1, se=TRUE, aes(group=fancy_strain, colour=fancy_strain), method = 'auto')+
  theme_bw(24)+
  #geom_jitter()+
  ylab("q90 EXT")+
  scale_x_discrete(labels=c("DMSO" = "DMSO", "6_25" = "6.25µM","12_5" = "12.5µM","25"="25µM","50"="50µM","100"="100µM"))+
  scale_fill_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y","1138"="A185P","1140"="A185P","1138"="Q131L","1098"="S145F", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V","1137"="Q131L","1139"="A185P","1314"="M257I","1317"="M257I","1097"="S145F"), values = cols)+
  scale_color_manual(name = "Strain", labels = c("N2" = "Suceptible","882"="Deletion","919"="F200Y","920"="F200Y","1138"="A185P","1140"="A185P","1138"="Q131L","1098"="S145F", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V","1137"="Q131L","1139"="A185P","1314"="M257I","1317"="M257I","1097"="S145F"), values = cols)+
  base_theme+
  theme(axis.title.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size=28,face = "bold"),
        axis.text.y=element_text(size=28,face = "bold"),
        axis.title.y=element_text(size=30,face = "bold"),
        legend.key.size = unit(0.4,"in"),
        legend.title = element_text(size=28,face = "bold"),
        legend.text = element_text(size=28,face = "bold"),
        legend.position = "None")
ggsave("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/plots/20190618_FBZdose/subtract_dose_q90ext.png", device = "png", plot = last_plot(), width = 14, height = 6, units = "in")

