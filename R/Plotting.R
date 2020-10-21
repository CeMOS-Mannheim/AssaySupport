#' get max value of df subsetted by x
#'
#' @param labels_df labels 
#' @param x subseting value
#'
#' @return max value
get_max <- function(labels_df, x) {
  max(labels_df[grepl(x, labels_df$species),"int"])
}


#' Overview of different conditions with annotations
#'
#' @param spec_ID     Integer, ID of first spectrum
#' @param n_overlay   Integer, number of spectra to plot (counting starts at \code{ID}) 
#' @param labels_df   Data.frame, df with labels
#' @param specs       \code{MALDIquant::MassSpectrum}
#' @param ab_max      Numeric, max intensity value for abeta region
#' @param aicd_max    Numeric, max intensity value for aicd region
#' @param title       Character, title of plot. If set to "auto" it is generated automatically
#' @param total_xlim  Numeric vector, x-axis limits for total spectrum view
#' @param abeta_xlim  Numeric vector, x-axis limits for abeta spectrum view
#' @param aicd_xlim   Numeric vector, x-axis limits for AICD spectrum view
#' @param substr_xlim Numeric vector, x-axis limits for substrate spectrum view
#' 
#' @importFrom magrittr %>%
#'
#' @return ggplot
#' @export
ggoverview <- function(spec_ID,
                       n_overlay,
                       labels_df,
                       specs = spectra_align, 
                       ab_max = get_max(labels_df, "Ab"),
                       aicd_max = get_max(labels_df, "AICD"),
                       title = "auto",
                       total_xlim = c(3000,16000),
                       abeta_xlim = c(4000,5000),
                       aicd_xlim = c(9400,10000),
                       substr_xlim = c(14500, 15500)) {
  
  
  total_spec <- ggoverlay_spectra(spec_ID = spec_ID, 
                                  specs = specs,
                                  labels_df = labels_df,
                                  xlim = total_xlim, 
                                  n_overlay = n_overlay, 
                                  normalizing = NA) 
  
  abeta_spec <- ggoverlay_spectra(spec_ID = spec_ID, 
                                  specs = specs,
                                  labels_df = labels_df,
                                  xlim = abeta_xlim, 
                                  n_overlay = n_overlay, 
                                  normalizing = NA) 
  if(!is.na(ab_max) && !ab_max == -Inf) {
    abeta_spec <- abeta_spec + scale_y_continuous(limits = c(0, ab_max))
  }
  
  aicd_spec <- ggoverlay_spectra(spec_ID = spec_ID, 
                                 specs = specs,
                                 labels_df = labels_df,
                                 xlim = aicd_xlim, 
                                 n_overlay = n_overlay, 
                                 normalizing = NA) + theme(legend.position = "none") + scale_x_continuous(breaks = seq(aicd_xlim[1], aicd_xlim[2], 250))
  if(!is.na(aicd_max)&& !aicd_max == -Inf) {
    aicd_spec <- aicd_spec + scale_y_continuous(limits = c(0, aicd_max))
  }
  
  substr_spec <- ggoverlay_spectra(spec_ID = spec_ID, 
                                   specs = specs,
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
  
  if(is.na(title)) {
    return(p)
  } else if(title == "auto") {
  p <- ggpubr::annotate_figure(p,
                               fig.lab = paste0(
                                 dplyr::filter(labels_df, plotIdx == spec_ID)$Substrate[1], 
                                 " incubated with ", 
                                 dplyr::filter(labels_df, plotIdx == spec_ID)$Type[1]))  
  } else {
    p <- ggpubr::annotate_figure(p,
                                 fig.lab = title)  
  }
  
  return(p)
}


#' Plot spectrum with annotations
#'
#' @param spec_ID       Integer, ID spectrum to plot
#' @param labels_df     Data.frame, df with labels
#' @param specs         \code{MALDIquant::MassSpectrum} 
#' @param xlim          Numeric vector, x-axis limits
#' @param n_overlay     Integer, number of spectra to plot (counting starts at \code{ID})
#' @param showLabel     Logical, show annotation in spectra
#' @param normalizing   Character or numeric, normalization to be applied to data. Useful to compare different plots.
#'                      Options:
#'                        "NA"     No normalization
#'                        "max"    Set highest signal to 100%
#'                        numeric  Normalize to peak at \code{normalizing +- tol}  
#' @param tol           Numeric, tolerance for normalization (if \code{normalizing} is numeric) 
#' @param separateInto  Character vector, splitt spectra names into columns. Sep = "_".
#' @importFrom magrittr %>% 
#' @importFrom tidyr separate
#' @importFrom ggplot2 aes ggplot geom_line 
#' @importFrom dplyr group_by
#'
#' @return ggplot
#' @export
#' 
#' @examples 
#' ggoverlay_spectra(spec_ID = 1, 
#'                   labels_df = test_result_df, 
#'                   specs = test_spectra_proc, 
#'                   n_overlay = 2)
ggoverlay_spectra <- function(spec_ID,
                              labels_df,
                              specs,
                              xlim = c(4000, 5000),
                              n_overlay = 3,
                              showLabel = T,
                              normalizing = NA,
                              tol = 3,
                              separateInto = c("Substrate", "Type", "Condition")) {
  if(is.null(names(specs))) {
    warning("Supplied spectra had no name. Setting name to noName.\n")
    names(specs) <- "noName"
  }
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
                    !mz > xlim[2]) 
  if(!any(is.na(separateInto))) {
    spec_df <- spec_df  %>% 
      separate(Sample, into = separateInto, sep = "_")
    separateInto <- "Sample"
  }
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
        left_join(spec_df_norm, by = separateInto) %>%
        group_by(Substrate, Type, Condition) %>%
        mutate(int = int/norm * 100)
      
      labels_df_norm <- labels_df %>%
        filter(!mz > (normalizing + tol),  !mz < (normalizing - tol)) %>%
        group_by(Substrate, Type, Condition) %>%
        summarise(norm = max(int)) %>%
        select(Substrate, Type, Condition, norm)
      
      labels_df <- labels_df %>%
        left_join(spec_df_norm, by = separateInto) %>%
        group_by(Substrate, Type, Condition) %>%
        mutate(int = int/norm * 100)
    }
  }
  
  
  
  
  p <- ggplot2::ggplot(spec_df, aes(x = mz, y = int, col = Condition)) + 
    ggplot2::geom_line(alpha = 0.8) +
    xtheme() 
  
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