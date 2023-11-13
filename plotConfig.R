theme_set(theme_bw()) #ggplot global setting
user='TurnerLinux'
#user='TurnerMac'
if (user == 'TurnerMac'){
	overleaf_plot_directory = '/Users/turnersilverthorn/research/overleaf/samplingPaper/figures/'
} else {
  overleaf_plot_directory = '/home/turner/research/overleaf/rate_limited_sampling/figures/'
}

control_labels <- function(plot, 
                           x_label = T, y_label = T, 
                           x_ticks = T, y_ticks = T,
                           legend_vis = T) {
  if (!x_label) {
    plot <- plot + theme(axis.title.x = element_blank())
  }
  
  if (!y_label) {
    plot <- plot + theme(axis.title.y = element_blank())
  }
  
  if (!x_ticks) {
    plot <- plot + theme(axis.text.x = element_blank())
  }
  
  if (!y_ticks) {
    plot <- plot + theme(axis.text.y = element_blank())
  }
  
  #TODO legend_vis option
  return(plot)
}