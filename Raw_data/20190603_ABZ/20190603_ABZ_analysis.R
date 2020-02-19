library(easysorter)
library(COPASutils)
library(dplyr)
library(tidyverse)
library(cowplot)

###Read in day 1 and analyze
dir <- c("~/Desktop/benzimidazoles/Albendazole/20190530_ABZV2/20190603_ABZ/")

raw <- read_data(dir)

raw_noncontam <- remove_contamination(raw)

raw_nongfp <- raw_noncontam[[1]]%>%
  mutate(row.col = paste(row,col))

mylist <-matrix(0,ncol=27,nrow = 0)

mylist <-data.frame(mylist)

colnames(mylist) <- c("date","experiment","round","assay","plate","row","col","sort","TOF","EXT","time","green","yellow","red","norm.EXT","norm.green","norm.yellow","norm.red","object","call50","stage","strain","condition","control","contamination","row.col","fit$classification")

summedplate <- sumplate(raw_noncontam, directories= FALSE, quantiles = TRUE)


biopruned <- bioprune(summedplate)

out_pruned <- prune_outliers(biopruned)

out_pruned_1 <- out_pruned %>%
 filter(strain != "1139")%>%
  mutate(assay=1)
save(out_pruned, file="non_regressed_ABZV2.RData")

out_parasite <- out_pruned_1 %>%
  filter(strain %in% parasite)
out_elegans <- out_pruned_1 %>%
  filter(strain %in% elegans)
out_figure <- out_pruned %>%
  filter(strain %in% parasite | strain %in% elegans)

full <- regress(out_pruned)
full_para <- regress(out_parasite)
full_ele <- regress(out_elegans)
full_1 <- regress(out_pruned_1)
full_fig <- regress(out_figure)

full_alleles <- full_1 %>%
  mutate(alleles = case_when(strain == "N2" ~ "WT",
                             strain == "882" ~ "Del",
                             strain == c("919", "920") ~ "F200Y", 
                             strain == c("1081", "1082") ~ "E198A", 
                             strain == "1327" ~ "E198V", 
                             strain == "1325" ~ "E198L",
                             strain == c("1075", "1076") ~ "F167Y",
                             strain == c("1137" , "1138") ~ "Q131L",
                             strain == c("1139" , "1140") ~ "A185P",
                             strain == "1317" ~ "M257I",
                             strain == c("1097" , "1098") ~ "S145F"))
cols <- c("N2" = "orange", "920" = "red", "1081" = "blue", "1075" = "green", "1325" = "brown","1138"="pink", "1098" = "maroon", "1327" = "yellow","N2" = "orange", "919" = "red", "1082" = "blue", "1076" = "green","1075" = "green", "1326" = "brown", "1328" = "yellow","1139"="black","1140"="black","1137"="pink","1317"="gray","1097"="maroon")

elegans <- c("N2", "882", "1137", "1140", "1097","1317")
parasite <- c("N2","882","1325","1327","1082","919","1076")
out_parasite %>%
  filter(trait=="mean.EXT",
         condition == "DMSO",
        # strain %in% parasite,
         !is.na(strain))%>%
         #strain %in% elegans) %>%
         #condition=="DMSO") %>%
         #strain==c("N2","882","1325","1327","1082","1097","919")) %>%
  dplyr::mutate(fancy_strain=factor(strain, levels = c("N2", "882","919","920","1081","1082","1327", "1325","1075","1076","1137","1138","1139","1140","1317","1097","1098")))%>%
  ggplot()+
  aes(x=fancy_strain, y = phenotype)+
  geom_boxplot(aes(fill=fancy_strain),outlier.shape = NA)+
  geom_jitter(width = 0.1)+
  #facet_grid(cols = vars(alleles),scales = "free_x")+
  theme_bw(24)+
 # ylim(-100,50)+
  ylab("Optical Density")+
  #scale_x_discrete(labels=c("N2" = "N2", "882" = "Del","919" = "F200Y","920"="F200Y","1139"="A185P","1081"="E198A","1138"="Q131L","1076"="F167Y", "1325" = "E198L", "1327" = "E198V", "1082" = "E198A","1075"="F167Y","1137"="Q131L","1097"="S145F","1098"="S145F","1140" = "A185P","1317"="M257I"))+
  scale_fill_manual(name = "Strain", labels = c("N2" = "N2","882"="Del","919"="F200Y","1140"="A185P","1138"="Q131L","1098"="S145F", "920" = "F200Y", "1081" = "E198A","1082" = "E198A", "1075" = "F167Y","1076" = "F167Y", "1325" = "E198L","1326" = "E198L", "1327" = "E198V","1328" = "E198V","1137"="Q131L","1314"="Q131L","1317"="M257I","1097"="S145F"), values = cols)+
  theme(axis.title.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size=20),
        axis.text.y=element_text(size=0),
        axis.title.y=element_text(size=30),
        legend.position = "None")
ggsave("full_DMSOmeanEXT_abz.pdf", device = "pdf", plot = last_plot(), width = 13, height = 8, units = "in") 
