load("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Processed_data/20190603regressed.Rdata" )

cols <- c("N2" = "orange", "920" = "red", "1081" = "blue", "1075" = "green", "1325" = "brown","1138"="pink", "1098" = "maroon", "1327" = "yellow","N2" = "orange", "919" = "red", "1082" = "blue", "1076" = "green","1075" = "green", "1326" = "brown", "1328" = "yellow","1139"="black","1140"="black","1137"="pink","1317"="gray","1097"="maroon")

full %>%
  filter(trait=="median.EXT",
         !is.na(strain))%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "882","919","920","1081","1082","1327", "1325","1075","1076","1137","1138","1139","1140","1317","1097","1098")))%>%
  ggplot()+
  aes(x=fancy_strain, y = phenotype)+
  geom_boxplot(aes(fill=fancy_strain),outlier.shape = NA)+
  geom_jitter(width = 0.1)+
  #facet_grid(cols = vars(alleles),scales = "free_x")+
  theme_bw(24)+
  # ylim(-100,50)+
  ylab("Median EXT")+
  scale_x_discrete(labels=c("N2" = "N2", "882" = "Del","919" = "F200Y","920"="F200Y","1139"="A185P","1081"="E198A","1138"="Q131L","1076"="F167Y", "1325" = "E198L", "1327" = "E198V", "1082" = "E198A","1075"="F167Y","1137"="Q131L","1097"="S145F","1098"="S145F","1140" = "A185P","1317"="M257I"))+
  scale_fill_manual(name = "Strain", labels = c("N2" = "N2","882"="Del","919"="F200Y","1140"="A185P","1138"="Q131L","1098"="S145F", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V","1137"="Q131L","1314"="Q131L","1317"="M257I","1097"="S145F"), values = cols)+
  theme(axis.title.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size=20),
        axis.text.y=element_text(size=0),
        axis.title.y=element_text(size=30),
        legend.position = "None")
ggsave("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/plots/20190603_ABZv2/regressed_medEXT.png", device = "png", plot = last_plot(), width = 14, height = 6, units = "in")

full %>%
  filter(trait=="median.TOF",
         !is.na(strain))%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "882","919","920","1081","1082","1327", "1325","1075","1076","1137","1138","1139","1140","1317","1097","1098")))%>%
  ggplot()+
  aes(x=fancy_strain, y = phenotype)+
  geom_boxplot(aes(fill=fancy_strain),outlier.shape = NA)+
  geom_jitter(width = 0.1)+
  #facet_grid(cols = vars(alleles),scales = "free_x")+
  theme_bw(24)+
  # ylim(-100,50)+
  ylab("Median TOF")+
  scale_x_discrete(labels=c("N2" = "N2", "882" = "Del","919" = "F200Y","920"="F200Y","1139"="A185P","1081"="E198A","1138"="Q131L","1076"="F167Y", "1325" = "E198L", "1327" = "E198V", "1082" = "E198A","1075"="F167Y","1137"="Q131L","1097"="S145F","1098"="S145F","1140" = "A185P","1317"="M257I"))+
  scale_fill_manual(name = "Strain", labels = c("N2" = "N2","882"="Del","919"="F200Y","1140"="A185P","1138"="Q131L","1098"="S145F", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V","1137"="Q131L","1314"="Q131L","1317"="M257I","1097"="S145F"), values = cols)+
  theme(axis.title.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size=20),
        axis.text.y=element_text(size=0),
        axis.title.y=element_text(size=30),
        legend.position = "None")
ggsave("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/plots/20190603_ABZv2/regressed_medTOF.png", device = "png", plot = last_plot(), width = 14, height = 6, units = "in")

full %>%
  filter(trait=="q90.TOF",
         !is.na(strain))%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "882","919","920","1081","1082","1327", "1325","1075","1076","1137","1138","1139","1140","1317","1097","1098")))%>%
  ggplot()+
  aes(x=fancy_strain, y = phenotype)+
  geom_boxplot(aes(fill=fancy_strain),outlier.shape = NA)+
  geom_jitter(width = 0.1)+
  #facet_grid(cols = vars(alleles),scales = "free_x")+
  theme_bw(24)+
  # ylim(-100,50)+
  ylab("q90 TOF")+
  scale_x_discrete(labels=c("N2" = "N2", "882" = "Del","919" = "F200Y","920"="F200Y","1139"="A185P","1081"="E198A","1138"="Q131L","1076"="F167Y", "1325" = "E198L", "1327" = "E198V", "1082" = "E198A","1075"="F167Y","1137"="Q131L","1097"="S145F","1098"="S145F","1140" = "A185P","1317"="M257I"))+
  scale_fill_manual(name = "Strain", labels = c("N2" = "N2","882"="Del","919"="F200Y","1140"="A185P","1138"="Q131L","1098"="S145F", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V","1137"="Q131L","1314"="Q131L","1317"="M257I","1097"="S145F"), values = cols)+
  theme(axis.title.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size=20),
        axis.text.y=element_text(size=0),
        axis.title.y=element_text(size=30),
        legend.position = "None")
ggsave("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/plots/20190603_ABZv2/regressed_q90TOF.png", device = "png", plot = last_plot(), width = 14, height = 6, units = "in")

full %>%
  filter(trait=="q90.EXT",
         !is.na(strain))%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "882","919","920","1081","1082","1327", "1325","1075","1076","1137","1138","1139","1140","1317","1097","1098")))%>%
  ggplot()+
  aes(x=fancy_strain, y = phenotype)+
  geom_boxplot(aes(fill=fancy_strain),outlier.shape = NA)+
  geom_jitter(width = 0.1)+
  #facet_grid(cols = vars(alleles),scales = "free_x")+
  theme_bw(24)+
  # ylim(-100,50)+
  ylab("q90 EXT")+
  scale_x_discrete(labels=c("N2" = "N2", "882" = "Del","919" = "F200Y","920"="F200Y","1139"="A185P","1081"="E198A","1138"="Q131L","1076"="F167Y", "1325" = "E198L", "1327" = "E198V", "1082" = "E198A","1075"="F167Y","1137"="Q131L","1097"="S145F","1098"="S145F","1140" = "A185P","1317"="M257I"))+
  scale_fill_manual(name = "Strain", labels = c("N2" = "N2","882"="Del","919"="F200Y","1140"="A185P","1138"="Q131L","1098"="S145F", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V","1137"="Q131L","1314"="Q131L","1317"="M257I","1097"="S145F"), values = cols)+
  theme(axis.title.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size=20),
        axis.text.y=element_text(size=0),
        axis.title.y=element_text(size=30),
        legend.position = "None")
ggsave("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/plots/20190603_ABZv2/regressed_q90EXT.png", device = "png", plot = last_plot(), width = 14, height = 6, units = "in")

full %>%
  filter(trait=="mean.TOF",
         !is.na(strain))%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "882","919","920","1081","1082","1327", "1325","1075","1076","1137","1138","1139","1140","1317","1097","1098")))%>%
  ggplot()+
  aes(x=fancy_strain, y = phenotype)+
  geom_boxplot(aes(fill=fancy_strain),outlier.shape = NA)+
  geom_jitter(width = 0.1)+
  #facet_grid(cols = vars(alleles),scales = "free_x")+
  theme_bw(24)+
  # ylim(-100,50)+
  ylab("mean TOF")+
  scale_x_discrete(labels=c("N2" = "N2", "882" = "Del","919" = "F200Y","920"="F200Y","1139"="A185P","1081"="E198A","1138"="Q131L","1076"="F167Y", "1325" = "E198L", "1327" = "E198V", "1082" = "E198A","1075"="F167Y","1137"="Q131L","1097"="S145F","1098"="S145F","1140" = "A185P","1317"="M257I"))+
  scale_fill_manual(name = "Strain", labels = c("N2" = "N2","882"="Del","919"="F200Y","1140"="A185P","1138"="Q131L","1098"="S145F", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V","1137"="Q131L","1314"="Q131L","1317"="M257I","1097"="S145F"), values = cols)+
  theme(axis.title.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size=20),
        axis.text.y=element_text(size=0),
        axis.title.y=element_text(size=30),
        legend.position = "None")
ggsave("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/plots/20190603_ABZv2/regressed_meanTOF.png", device = "png", plot = last_plot(), width = 14, height = 6, units = "in")

full %>%
  filter(trait=="mean.EXT",
         !is.na(strain))%>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "882","919","920","1081","1082","1327", "1325","1075","1076","1137","1138","1139","1140","1317","1097","1098")))%>%
  ggplot()+
  aes(x=fancy_strain, y = phenotype)+
  geom_boxplot(aes(fill=fancy_strain),outlier.shape = NA)+
  geom_jitter(width = 0.1)+
  #facet_grid(cols = vars(alleles),scales = "free_x")+
  theme_bw(24)+
  # ylim(-100,50)+
  ylab("mean EXT")+
  scale_x_discrete(labels=c("N2" = "N2", "882" = "Del","919" = "F200Y","920"="F200Y","1139"="A185P","1081"="E198A","1138"="Q131L","1076"="F167Y", "1325" = "E198L", "1327" = "E198V", "1082" = "E198A","1075"="F167Y","1137"="Q131L","1097"="S145F","1098"="S145F","1140" = "A185P","1317"="M257I"))+
  scale_fill_manual(name = "Strain", labels = c("N2" = "N2","882"="Del","919"="F200Y","1140"="A185P","1138"="Q131L","1098"="S145F", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V","1137"="Q131L","1314"="Q131L","1317"="M257I","1097"="S145F"), values = cols)+
  theme(axis.title.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size=20),
        axis.text.y=element_text(size=0),
        axis.title.y=element_text(size=30),
        legend.position = "None")
ggsave("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/plots/20190603_ABZv2/regressed_meanEXT.png", device = "png", plot = last_plot(), width = 14, height = 6, units = "in")

