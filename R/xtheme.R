#' Standard ggplot theme
#' @export
xtheme <- function() 
{
  ggplot2::theme(legend.position = "right", legend.title = element_text(face = "bold", 
                                                                        colour = "#000000", size = 10), legend.background = element_rect(fill = "#FFFFFF"), 
                 legend.key = element_rect(fill = "#FFFFFF", colour = "#FFFFFF"), 
                 legend.text = element_text(colour = "#000000", size = 10), 
                 plot.title = element_text(colour = "black", face = "bold", 
                                           size = 18, vjust = 1), plot.background = element_rect(fill = "white", 
                                                                                                 colour = "white"), panel.background = element_rect(fill = "white"), 
                 axis.text = element_text(colour = "#000000", size = 14), 
                 axis.text.x = element_text(colour = "#000000", size = 12), 
                 axis.title = element_text(colour = "black", face = "bold", 
                                           size = 12), axis.title.y = element_text(colour = "black", 
                                                                                   face = "bold", size = 12, vjust = 1), axis.ticks = element_line(colour = "black"), 
                 panel.grid.major.x = element_line(colour = "#FFFFFF"), 
                 panel.grid.minor.x = element_line(colour = "#FFFFFF"), 
                 panel.grid.major.y = element_line(colour = "#FFFFFF"), 
                 panel.grid.minor.y = element_line(colour = "#FFFFFF"), 
                 strip.text = element_text(colour = "white"), strip.background = element_rect(fill = "#333333"))
}
