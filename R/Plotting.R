#' get max value of df subsetted by x
#'
#' @param labels_df labels 
#' @param x subseting value
#'
#' @return max value
get_max <- function(labels_df,
                    x) {
  max(labels_df[grepl(x, labels_df$species),"int"])
}


#' Overview using ggplot
#'
#' @param spec_ID ID
#' @param n_overlay number of spectra 
#' @param labels_df df with labels
#' @param specs spectra
#' @param ab_max max value for abeta region
#' @param aicd_max max value for aicd region
#' @importFrom magrittr %>%
#'
#' @return ggplot
#' @export
ggoverview <- function(spec_ID,
                       n_overlay,
                       labels_df = label_df,
                       specs = spectra_align, 
                       ab_max = get_max(labels_df, "Ab"),
                       aicd_max = get_max(labels_df, "AICD"),
                       label = "auto",
                       total_xlim = c(3000,16000),
                       abeta_xlim = c(4000,5000),
                       aicd_xlim = c(9400,10000),
                       substr_xlim = c(14500, 15500)) {
  
  
  total_spec <- ggoverlay_spectra(spec_ID = spec_ID, 
                                  labels_df = labels_df,
                                  xlim = total_xlim, 
                                  n_overlay = n_overlay, 
                                  normalizing = NA) 
  
  abeta_spec <- ggoverlay_spectra(spec_ID = spec_ID, 
                                  labels_df = labels_df,
                                  xlim = abeta_xlim, 
                                  n_overlay = n_overlay, 
                                  normalizing = NA) 
  if(!is.na(ab_max) && !ab_max == -Inf) {
    abeta_spec <- abeta_spec + scale_y_continuous(limits = c(0, ab_max))
  }
  
  aicd_spec <- ggoverlay_spectra(spec_ID = spec_ID, 
                                 labels_df = labels_df,
                                 xlim = aicd_xlim, 
                                 n_overlay = n_overlay, 
                                 normalizing = NA) + theme(legend.position = "none") + scale_x_continuous(breaks = seq(aicd_xlim[1], aicd_xlim[2], 250))
  if(!is.na(aicd_max)&& !aicd_max == -Inf) {
    aicd_spec <- aicd_spec + scale_y_continuous(limits = c(0, aicd_max))
  }
  
  substr_spec <- ggoverlay_spectra(spec_ID = spec_ID, 
                                   labels_df = labels_df,
                                   xlim = substr_xlim, 
                                   n_overlay = n_overlay, 
                                   normalizing = NA) + theme(legend.position = "none") + scale_x_continuous(breaks = seq(substr_xlim[1], substr_xlim[2], 500))
  
  p <- ggpubr::ggarrange(total_spec, 
                         abeta_spec, 
                         ggpubr::ggarrange(aicd_spec, 
                                           substr_spec, 
                                           nrow = 1, ncol = 2), 
                         common.legend = TRUE, 
                         legend = "bottom",
                         ncol = 1, nrow = 3)
  
  if(is.na(label)) {
    return(p)
  } else if(label == "auto") {
  p <- ggpubr::annotate_figure(p,
                               fig.lab = paste0(
                                 dplyr::filter(labels_df, plotIdx == spec_ID)$Substrate[1], 
                                 " incubated with ", 
                                 dplyr::filter(labels_df, plotIdx == spec_ID)$Type[1]))  
  } else {
    p <- ggpubr::annotate_figure(p,
                                 fig.lab = label)  
  }
  
  return(p)
}


#' Plot spectrum with ggplot
#'
#' @param spec_ID ID
#' @param labels_df df with labels
#' @param specs spectra 
#' @param xlim limits
#' @param n_overlay n
#' @param showLabel logical
#' @param normalizing logical
#' @param tol tolerance
#' @importFrom magrittr %>%
#'
#' @return ggplot
#' @export
ggoverlay_spectra <- function(spec_ID,
                              labels_df = label_df,
                              specs = spectra_align,
                              xlim = c(4000, 5000),
                              n_overlay = 3,
                              showLabel = T,
                              normalizing = "max",
                              tol = 3) {
  spec_df <- data.frame(Sample = character(),
                        mz = numeric(),
                        int = numeric())
  
  for(i in 0:(n_overlay-1)) {
    df <- data.frame(Sample = names(specs[spec_ID+i]),
                     mz = specs[[spec_ID+i]]@mass,
                     int = specs[[spec_ID+i]]@intensity)
    spec_df <- rbind(spec_df, df)
  }
  
  spec_df <- filter(spec_df, 
                    !mz < xlim[1], 
                    !mz > xlim[2]) %>%
    separate(Sample, into = c("Substrate", "Type", "Condition"), sep = "_")
  
  if(!is.na(normalizing)) {
    if(normalizing == "max") {
      spec_df <- spec_df %>%
        group_by(Substrate, Type, Condition) %>%
        mutate(int = int/max(int) * 100)
      
      labels_df <- labels_df %>%
        filter(!mz < xlim[1],
               !mz > xlim[2]) %>%
        group_by(Substrate, Type, Condition) %>%
        mutate(int = int/max(int) * 100)
    }
    if(is.numeric(normalizing)) {
      spec_df_norm <- spec_df %>%
        filter(!mz > (normalizing + tol),  !mz < (normalizing - tol)) %>%
        group_by(Substrate, Type, Condition) %>%
        summarise(norm = max(int)) %>%
        select(Substrate, Type, Condition, norm)
      
      spec_df <- spec_df %>%
        left_join(spec_df_norm, by =c("Substrate", "Type", "Condition")) %>%
        group_by(Substrate, Type, Condition) %>%
        mutate(int = int/norm * 100)
      
      labels_df_norm <- labels_df %>%
        filter(!mz > (normalizing + tol),  !mz < (normalizing - tol)) %>%
        group_by(Substrate, Type, Condition) %>%
        summarise(norm = max(int)) %>%
        select(Substrate, Type, Condition, norm)
      
      labels_df <- labels_df %>%
        left_join(spec_df_norm, by =c("Substrate", "Type", "Condition")) %>%
        group_by(Substrate, Type, Condition) %>%
        mutate(int = int/norm * 100)
    }
  }
  
  
  
  
  p <- ggplot2::ggplot(spec_df, aes(x = mz, y = int, col = Condition)) + 
    ggplot2::geom_line(alpha = 0.8) +
    AlzTools::xtheme() 
  
  suppressWarnings(
    if(!is.na(labels_df)) {
      labels_df <- filter(labels_df, 
                          !mz < xlim[1],
                          !mz > xlim[2],
                          Substrate %in% spec_df$Substrate,
                          Type %in% spec_df$Type,
                          Condition %in% spec_df$Condition) 
      
      p <- p + ggplot2::geom_point(data = labels_df, aes(x = mz, y = int), alpha = 0.8)
      
      if(showLabel) {
        l_df <- labels_df %>% 
          group_by(species) %>%
          filter(!is.na(species),
                 dplyr::row_number() == 1)
        p <- p + ggrepel::geom_text_repel(data = l_df, 
                                 aes(x = mz, y =  int, label = species), 
                                 col = "black", 
                                 min.segment.length = 1, 
                                 nudge_y = max(l_df$int)*0.15)
      }
    })
  
  if(!is.na(normalizing)) {
    
    if(normalizing == "max") {
      p <- p + labs(y = "Rel. Intensity as %\nof max(Int.)")
    } else if(is.numeric(normalizing)) {
      p <- p + labs(y = paste0("Rel. Intensity as %\nof Int. at mz ", normalizing))
    }
  }
  return(p)
}

#' Overview using base R graphics
#'
#' @param idx index
#' @param specs spectra
#' @param labeldf df with labels
#' @param total_xlim limits
#' @param AICD_xlim limits
#' @param Ab_xlim limits
#' @param subst_xlim limits
#' @param offset offset
#' @param ymax_label n
#' @param ymax_adjust n
#' @param overview_adjust n
#'
#' @return plot
#' @export
overviewplot <- function(idx, 
                         specs = spectra_align, 
                         labeldf = label_df, 
                         total_xlim = c(3000, 16000), 
                         AICD_xlim = c(9500, 10000),
                         Ab_xlim =  c(3500, 5500),
                         subst_xlim = c(13500, 15500),
                         offset = 1.9,
                         ymax_label = TRUE,
                         ymax_adjust = 1.7,
                         overview_adjust = 5) {
  l_df <- labeldf[labeldf$plotIdx == idx,] 
  if(ymax_label) {
    l_df <- filter(l_df, !is.na(species))
    if(dim(l_df)[1] < 1)
      l_df <- labeldf[labeldf$plotIdx == idx,]
  }
  layout(matrix(c(1,1,1,1,2,2,2,2,3,3,4,4), 3, 4, byrow = TRUE))
  par(mar = c(1,1,1,0) + 1.2)
  plot(spectra_align[[idx]], xlim= total_xlim, ylim = c(0,max(l_df$int)*ymax_adjust*overview_adjust), main = metaData(spectra_align[[idx]])$sampleName)
  points(peaks[[idx]], col = "red")
  text(x = l_df$mz, y = l_df$int, labels = l_df$species, pos = 3, offset = offset, srt = 90)
  plot(spectra_align[[idx]], xlim= Ab_xlim, main = "Abeta-region", ylim = c(0,max(l_df$int*ymax_adjust)))
  points(peaks[[idx]], col = "red")
  text(x = l_df$mz, y = l_df$int, labels = l_df$species, pos = 3, offset = offset, srt = 90)
  plot(spectra_align[[idx]], xlim= AICD_xlim, main = "AICD-region", ylim = c(0,max(l_df$int*0.4)))
  points(peaks[[idx]], col = "red")
  text(x = l_df$mz, y = l_df$int, labels = l_df$species, pos = 3, offset = offset, srt = 90)
  plot(spectra_align[[idx]], xlim= subst_xlim, main = "Substrate-region", ylim = c(0,max(l_df$int*0.1)))
  points(peaks[[idx]], col = "red")
  text(x = l_df$mz, y = l_df$int, labels = l_df$species, pos = 3, offset = offset, srt = 90)
}

#' Plot comparing different mutations using base R graphics
#'
#' @param idx n
#' @param plot_n n
#' @param specs n
#' @param labeldf n
#' @param Ab_xlim n
#' @param ymax_label n
#' @param ymax_adjust n
#'
#' @return plot
#' @export
mutation_plot <- function(idx, 
                          plot_n = 4, 
                          specs = spectra_align, 
                          labeldf = label_df, 
                          Ab_xlim =  c(3500, 5500), 
                          ymax_label = TRUE, 
                          ymax_adjust = 1.7) {
  layout(matrix(1:plot_n, plot_n, 1, byrow = TRUE))
  par(mar = c(1,1,1,0) + 1.2)
  l_df <- labeldf[labeldf$plotIdx == idx,]
  if(ymax_label) {
    l_df <- filter(l_df, !is.na(species))
    if(dim(l_df)[1] < 1)
      l_df <- labeldf[labeldf$plotIdx == idx+n,]
  }
  plot(spectra_align[[idx]], xlim= Ab_xlim, ylim = c(0,max(l_df$int)*ymax_adjust), main = metaData(spectra_align[[idx]])$sampleName)
  points(peaks[[idx]], col = "red")
  text(x = l_df$mz, y = l_df$int, labels = l_df$species, pos = 3, offset = 0.5)
  for(n in 1:(plot_n-1)) {
    l_df <- labeldf[labeldf$plotIdx == idx+n,]
    if(ymax_label) {
      l_df <- filter(l_df, !is.na(species))
      if(dim(l_df)[1] < 1)
        l_df <- labeldf[labeldf$plotIdx == idx+n,]
    }
    plot(spectra_align[[idx+n]], xlim= Ab_xlim, ylim = c(0,max(l_df$int)*ymax_adjust), main = metaData(spectra_align[[idx+n]])$sampleName)
    points(peaks[[idx+n]], col = "red")
    text(x = l_df$mz, y = l_df$int, labels = l_df$species, pos = 3, offset = 0.5)
  }
}