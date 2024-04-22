theme_set(theme_classic()) #ggplot global setting
#user='TurnerLinux'
user='TurnerMac'
if (user == 'TurnerMac'){
	overleaf_plot_directory = '/Users/turnersilverthorn/research/overleaf/samplingPaper/figures/'
} else {
  overleaf_plot_directory = '/home/turner/research/ms_powerCHORD/figures/'
}

get_legend <- function(p) {
   tmp <- ggplot_gtable(ggplot_build(p))
   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
   legend <- tmp$grobs[[leg]]
   legend
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

fs_glob=9

rad_brk = c(0,pi/2,pi,3*pi/2,2*pi)
rad_lab = c(expression(0),
            expression(pi/2),
            expression(pi),
            expression(3*pi/2),
            expression(2*pi))
#  theme(
#   strip.background=element_blank(),
#   text=element_text(size=9),
#   plot.margin=margin(0,0,0,0))
#   panel.grid.major = element_blank(),panel.grid.minor = element_blank()
#)
