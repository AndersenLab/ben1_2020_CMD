source("~/Dropbox/AndersenLab/LabFolders/Clay/2019_benzimidazole_paper/Scripts/theme.R")
library(cowplot)
prune_outliers_CDedit <- function(data, iqr = FALSE) {
  
  # Make sure that the data being fed into the pruning function is in long
  # format
  
  data <- ensure_long_CDedit(data)
  
  # remove NA phenotypes
  
  napheno <- data[is.na(data[, "phenotype"]), ]
  
  # if type is IQR: 
  if(iqr == TRUE) {
    
    datacuts <- data %>%
      
      # Filter out all of the wash and/or empty wells
      dplyr::filter(!is.na(strain)) %>%
      
      # Group on condition, trait, and strain
      dplyr::group_by(condition, trait, strain) %>%
      
      # calculate the median, mean, 2IQR and 2sd
      dplyr::mutate(med = median(phenotype), 
                    iqr = IQR(phenotype),
                    mph = mean(phenotype), 
                    sph = sd(phenotype)) %>%
      dplyr::ungroup() %>%
      
      # Filter phenotype +/- 2*IQR
      dplyr::filter(phenotype >= med - 2*iqr, phenotype <= 2*iqr + med) %>%
      
      # select out columns not needed
      dplyr::select(-c(med:sph))
    
    # if type is not IQR, it should be SD
  } else {
    
    datacuts <- data %>%
      
      # Filter out all of the wash and/or empty wells
      dplyr::filter(!is.na(strain), !is.na(phenotype)) %>%
      
      # Group on condition, trait, and strain
      dplyr::group_by(condition, trait, strain) %>%
      
      # calculate the median, mean, 2IQR and 2sd
      dplyr::mutate(med = median(phenotype), 
                    iqr = IQR(phenotype),
                    mph = mean(phenotype), 
                    sph = sd(phenotype)) %>%
      dplyr::ungroup() %>%
      
      # Filter phenotype +/- 2*sd
      dplyr::filter(phenotype >= mph - 2*sph, phenotype <= 2*sph + mph) %>%
      
      # select out columns not needed
      dplyr::select(-c(med:sph))
  }
  
  return(datacuts)
}

ensure_long_CDedit <- function(data){
  if("trait" %in% colnames(data)){
    return(data)
  } else {
    longdata <- tidyr::gather(data, trait, phenotype, -c(experiment,plate, condition, strain, row, col))
    return(longdata)
  }
}

plot_manhatten <- function(data, option = "top") {
  
  if (isTRUE(option == "top")){
  data %>%
    dplyr::filter(CHROM != "MtDNA") %>%
    dplyr::distinct(marker, .keep_all = T) %>%
    dplyr::mutate(EIGEN_CUTOFF = -log10(.05/772)) %>%
    dplyr::mutate(EIGEN_SIG = ifelse(log10p > BF, "1", 
                                     ifelse(log10p > EIGEN_CUTOFF, "2", "0")) )%>%
    ggplot2::ggplot(.) +
    ggplot2::aes(x = POS/1e6, y = log10p) +
    ggplot2::scale_color_manual(values = c("0" = "black", 
                                           "1" = "red",
                                           "2" = "hotpink3")) +
    ggplot2::scale_alpha_manual(values = c("0" = 0.5, 
                                           "1" = 1,
                                           "2" = 1)) +
    #ggplot2::geom_rect(ggplot2::aes(xmin = startPOS/1e6, 
     #                               xmax = endPOS/1e6, 
      #                              ymin = 0, 
       #                             ymax = Inf, 
        #                            fill = "blue"), 
         #              color = "blue",fill = "cyan",linetype = 2, 
          #             alpha=.3)+
    ggplot2::geom_hline(ggplot2::aes(yintercept = BF),
                        color = "gray60", 
                        alpha = .75,  
                        size = 1) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = EIGEN_CUTOFF),
                        color = "gray60", 
                        alpha = .75,  
                        size = 1,
                        linetype = 2) +
    ggplot2::geom_point( ggplot2::aes(color= factor(EIGEN_SIG), alpha = factor(EIGEN_SIG)) ) +
    ggplot2::facet_grid( . ~ CHROM, scales = "free_x" , space = "free_x") +
    ggplot2::theme_bw(20) +
    scale_y_continuous() +
    base_theme+
    ggplot2::theme(strip.background = element_blank(),
                   legend.position = "none",
                   panel.grid.minor = element_blank(),
                   panel.grid.major = element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.title.x = element_blank(),
                   axis.text.x = element_blank(),
                   axis.text.y = element_text(size=5.5),
                   plot.margin = unit(c(0.1, 0.1, 0, 0.1), "in")) +
    ggplot2::labs(x = "Genomic Position (Mb)",
                  y = expression(-log[10](italic(p))))
  }
  else if (isTRUE(option == "mid")){
    data %>%
      dplyr::filter(CHROM != "MtDNA") %>%
      dplyr::distinct(marker, .keep_all = T) %>%
      dplyr::mutate(EIGEN_CUTOFF = -log10(.05/772)) %>%
      dplyr::mutate(EIGEN_SIG = ifelse(log10p > BF, "1", 
                                       ifelse(log10p > EIGEN_CUTOFF, "2", "0")) )%>%
      ggplot2::ggplot(.) +
      ggplot2::aes(x = POS/1e6, y = log10p) +
      ggplot2::scale_color_manual(values = c("0" = "black", 
                                             "1" = "red",
                                             "2" = "hotpink3")) +
      ggplot2::scale_alpha_manual(values = c("0" = 0.5, 
                                             "1" = 1,
                                             "2" = 1)) +
      #ggplot2::geom_rect(ggplot2::aes(xmin = startPOS/1e6, 
       #                               xmax = endPOS/1e6, 
        #                              ymin = 0, 
         #                             ymax = Inf, 
          #                            fill = "blue"), 
           #              color = "blue",fill = "cyan",linetype = 2, 
            #             alpha=.3)+
      ggplot2::geom_hline(ggplot2::aes(yintercept = BF),
                          color = "gray60", 
                          alpha = .75,  
                          size = 1) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = EIGEN_CUTOFF),
                          color = "gray60", 
                          alpha = .75,  
                          size = 1,
                          linetype = 2) +
      ggplot2::geom_point( ggplot2::aes(color= factor(EIGEN_SIG), alpha = factor(EIGEN_SIG)) ) +
      ggplot2::facet_grid( . ~ CHROM, scales = "free_x" , space = "free_x") +
      ggplot2::theme_bw(20) +
      scale_y_continuous() +
      base_theme+
      ggplot2::theme(strip.background = element_blank(),
                     axis.title.x = element_blank(),
                     strip.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.text.x=element_blank(),
                     legend.position = "none",
                     panel.grid.minor = element_blank(),
                     panel.grid.major = element_blank(),
                     plot.margin = unit(c(0.1, 0.1, 0, 0.1), "in")) +
      ggplot2::labs(x = "Genomic Position (Mb)",
                    y = expression(-log[10](italic(p))))
  }
  else{
    data %>%
      dplyr::filter(CHROM != "MtDNA") %>%
      dplyr::distinct(marker, .keep_all = T) %>%
      dplyr::mutate(EIGEN_CUTOFF = -log10(.05/772)) %>%
      dplyr::mutate(EIGEN_SIG = ifelse(log10p > BF, "1", 
                                       ifelse(log10p > EIGEN_CUTOFF, "2", "0")) )%>%
      ggplot2::ggplot(.) +
      ggplot2::aes(x = POS/1e6, y = log10p) +
      ggplot2::scale_color_manual(values = c("0" = "black", 
                                             "1" = "red",
                                             "2" = "hotpink3")) +
      ggplot2::scale_alpha_manual(values = c("0" = 0.5, 
                                             "1" = 1,
                                             "2" = 1)) +
      #ggplot2::geom_rect(ggplot2::aes(xmin = startPOS/1e6, 
       #                               xmax = endPOS/1e6, 
        #                              ymin = 0, 
         #                             ymax = Inf, 
          #                            fill = "blue"), 
           #              color = "blue",fill = "cyan",linetype = 2, 
            #             alpha=.3)+
      ggplot2::geom_hline(ggplot2::aes(yintercept = BF),
                          color = "gray60", 
                          alpha = .75,  
                          size = 1) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = EIGEN_CUTOFF),
                          color = "gray60", 
                          alpha = .75,  
                          size = 1,
                          linetype = 2) +
      ggplot2::geom_point( ggplot2::aes(color= factor(EIGEN_SIG), alpha = factor(EIGEN_SIG)) ) +
      ggplot2::facet_grid( . ~ CHROM, scales = "free_x" , space = "free_x") +
      ggplot2::theme_bw(20) +
      scale_y_continuous() +
      base_theme+
      ggplot2::theme(strip.background = element_blank(),
                     strip.text.x = element_blank(),
                     legend.position = "none",
                     panel.grid.minor = element_blank(),
                     panel.grid.major = element_blank(),
                     plot.margin = unit(c(0.1, 0.1, 0, 0.1), "in")) +
      ggplot2::labs(x = "Genomic Position (Mb)",
                    y = expression(-log[10](italic(p))))
  }
}

