## Plotting function for signal drift
cplt = function(x, dat_train, testing){

  dat_lab = dat_train
  list2env(x ,.GlobalEnv)

  df_pseudoQC = df_pseudoQC %>%
    mutate(area = ifelse(area==0, NA, area))
  df_pseudoQC_corrected = df_pseudoQC_corrected %>%
    mutate(area = ifelse(area==0, NA, area))

  class_qc = df_pseudoQC %>%
    group_by(batch) %>%
    mutate(index = 1:n()) %>%
    filter(class=="QC") %>%
    # mutate(area = log2(area)) %>%
    drop_na(area)

  class_pseudoqc = df_pseudoQC %>%
    group_by(batch) %>%
    mutate(index = 1:n()) %>%
    filter(class=="Pseudo_QC") %>%
    # mutate(area = log2(area)) %>%
    drop_na(area)

  metab = unique(class_pseudoqc$compound)

  pseudo_plot = df_pseudoQC %>%
    # mutate(area = log2(area)) %>%
    group_by(batch) %>%
    mutate(index = 1:n()) %>%
    ggplot(., aes(x = index, y = area, col = batch, shape = class)) +
    geom_point(size = 0.5)+
    geom_point(data = class_qc, size = 1.75, color = "black") +
    geom_point(data = class_pseudoqc, size = 1.75, shape = 19, show.legend = FALSE, color = "red") +
    scale_shape_manual(values=c(19, 8, 19))+
    ggsci::scale_color_uchicago()+
    geom_line(data = class_qc, aes(x = index, y = area), color = "black")+
    geom_line(data = class_pseudoqc, aes(x = index, y = area), color = "red")+
    scale_x_continuous(expand=expansion(mult=c(0.01, 0.03)))+
    facet_grid(cols = vars(batch), scales = "free_x", space = "free")+
    theme_classic()+
    theme(
      axis.text.x = element_text(vjust = 0.5, hjust=1, angle = 90, size = 12,face = "bold", color = "black"),
      axis.text.y = element_text(size = 12,face = "bold", color = "black"),
      axis.title.y = element_text(size = 12, color = "black", face = "bold"),
      axis.title.x = element_text(size = 16, color = "black", face = "bold"),
      strip.text.x = element_text(size = 14, color = "black", face = "bold"),
      plot.title = element_text(hjust = 0.5,color = "black", face = "bold"),
      panel.spacing.x = unit(0.9, "lines"),
      legend.position="top",
      legend.box="vertical",
      legend.margin=margin(),
      legend.justification='center',
      legend.title = element_blank(),
      legend.direction='horizontal'
    )+
    guides(colour = guide_legend(nrow = 1,override.aes = list(size=5)),
           shape = guide_legend(nrow = 1,override.aes = list(size=c(5,4,2))))+
    labs(x = "Order in run", y = paste0("log2(", testing, "+1)"), title = paste0(metab," raw data using ",dat_lab, " data to train pseudo-QC"))

  # legend
  legend <- cowplot::get_legend(pseudo_plot)

  class_qc1 = df_pseudoQC_corrected %>%
    group_by(batch) %>%
    mutate(index = 1:n()) %>%
    filter(class=="QC")
    # mutate(area_corrected = log2(area_corrected))
  class_pseudoqc1 = df_pseudoQC_corrected %>%
    group_by(batch) %>%
    mutate(index = 1:n()) %>%
    filter(class=="Pseudo_QC")
    # mutate(area_corrected = log2(area_corrected))

  pseudo_plot1 = df_pseudoQC_corrected %>%
    group_by(batch) %>%
    mutate(index = 1:n()) %>%
    # mutate(area_corrected = log2(area_corrected)) %>%
    ggplot(., aes(x = index, y = area_corrected, col = batch, shape = class)) +
    geom_point(size = 0.5)+
    geom_point(data = class_qc1, size = 1.75, color = "black") +
    geom_point(data = class_pseudoqc1, size = 1.75, shape = 19, show.legend = FALSE, color = "red") +
    scale_shape_manual(values=c(19, 8, 19))+
    ggsci::scale_color_uchicago()+
    geom_line(data = class_qc1, aes(x = index, y = area_corrected), color = "black")+
    geom_line(data = class_pseudoqc1, aes(x = index, y = area_corrected), color = "red")+
    scale_x_continuous(expand=expansion(mult=c(0.01, 0.03)))+
    facet_grid(cols = vars(batch), scales = "free_x", space = "free")+
    theme_classic()+
    theme(
      axis.text.x = element_text(vjust = 0.5, hjust=1, angle = 90, size = 12,face = "bold", color = "black"),
      axis.text.y = element_text(size = 12,face = "bold", color = "black"),
      axis.title.y = element_text(size = 12, color = "black", face = "bold"),
      axis.title.x = element_text(size = 16, color = "black", face = "bold"),
      strip.text.x = element_text(size = 14, color = "black", face = "bold"),
      plot.title = element_text(hjust = 0.5,color = "black", face = "bold"),
      panel.spacing.x = unit(0.9, "lines"),
      legend.position="top",
      legend.box="vertical",
      legend.margin=margin(),
      legend.justification='center',
      legend.title = element_blank(),
      legend.direction='horizontal'
    )+
    guides(colour = guide_legend(nrow = 2,override.aes = list(size=5)),
           shape = guide_legend(nrow = 1,override.aes = list(size=c(5,4,2))))+
    labs(x = "Order in run",y = paste0("log2(", testing, "+1)"), title = paste0(metab," signal drift corrected data"))

  c_plot = cowplot::plot_grid(legend, rel_heights = c(0.2,1,1),
                              pseudo_plot + theme(legend.position = "none"),
                              pseudo_plot1 + theme(legend.position = "none"),
                              ncol = 1, labels = c("","A","B"))


  return(c_plot)
}
