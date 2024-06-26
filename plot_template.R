require(ggplot2)
require(ggplotify)
require(patchwork)
require(dplyr)

# simulate data
N=1e2
rowvals = c('row1','row2')
colvals = c('col1','col2','col3')
df =data.frame(xvar=2*pi*runif(N),yvar=runif(N),zvar=runif(N),
               rvar=sample(rowvals,N,replace = T),
               cvar=sample(colvals,N,replace = T))

# rendered plot size
plt_height = 3
plt_width  = 6

# plot settings
fsize       = 10
xmin        = 0
xmax        = 2*pi
ymin        = 0
ymax        = 1
legend_vert = T
mute_axes   = F
radial_x    = F

# convert to factors
df$rvar = factor(df$rvar,c('row1','row2'),c('dog','cat'))
df$cvar = factor(df$cvar,c('col1','col2','col3'),c('temp1','temp2','temp3'))

# raw plot command
plt = df %>%  ggplot(aes(x=xvar,y=yvar,color=zvar))+
  geom_point()+
  facet_grid(rvar~cvar)

# theme
theme_set(theme_classic()) 

# font size
plt=plt+theme(text=element_text(size=fsize))


# axis labels 
plt = plt + labs(x=element_text('x variable'),
                 y=element_text('y variable'),
                 color='z variable')
if (mute_axes){
  plt = plt + theme(axis.title.x=element_blank()) 
  plt = plt + theme(axis.title.y=element_blank()) 
  plt = plt + theme(axis.text.x=element_blank()) 
  plt = plt + theme(axis.text.y=element_blank()) 
}

# axis ticks/limits (option to just do start/end start/end/mid)
plt = plt + scale_x_continuous(limits =c(xmin,xmax),
                               breaks=c(xmin,xmax),
                               labels=round(c(xmin,xmax),2)) 
plt = plt + scale_y_continuous(limits =c(ymin,ymax),
                               breaks=c(ymin,ymax),
                               labels=round(c(ymin,ymax),2)) 
if (radial_x){
  rad_brk = c(0,pi/2,pi,3*pi/2,2*pi)
  rad_lab = c(expression(0),expression(pi/2),
              expression(pi),expression(3*pi/2),
              expression(2*pi))
  plt = plt + scale_x_continuous(limits=c(0,2*pi),
                                 breaks =rad_brk,
                                 labels = rad_lab)
}

# color bar map/size/tweaks
plt=plt+scale_color_gradient2(
                              limits   = c(0,1),
                              midpoint = 0.5,
                              low=rgb(0.36, 0.54, 0.66),
                              high=rgb(.81,.1,.26),
                              )

if (legend_vert){
  plt=plt+guides(fill=guide_colorbar(title.position='right'))
  plt=plt+theme(legend.key.height = unit(plt_height*.15, "in"),
                 legend.title = element_text(angle = 90),
                 legend.title.align = 0.5,
                 legend.direction = "vertical")
}else{
  plt=plt+guides(fill=guide_colorbar(title.position='top'))
  plt=plt+theme(legend.position='bottom',
                 legend.key.width = unit(plt_width*.15, "in"),
                 legend.title.align = 0.5,
                 legend.direction = "horizontal")
}

# misc visuals
plt = plt + theme(
  strip.background=element_blank(),
  plot.margin = margin(0,0,0,0),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.text.x = element_text(vjust = 0.25)
  )

show_temp_plt=function(plt,plt_width,plt_height){
  plt_path <- tempfile(fileext = ".png")
  ggsave(plt_path, plt, width =plt_width, height = plt_height, units = "in",
         dpi = 96)
  
  viewer <- getOption("viewer")
  viewer(plt_path)
}

# how to get nice labels for wrapped variables
df=data.frame(val=runif(10),Nval=sample(c('1','2'),10,T)) 
df$Nval = as.factor(df$Nval)
df %>% 
  mutate(N=Nval) %>% 
  ggplot(aes(x=val,y=1))+geom_point()+
  facet_wrap(~N,labeller = purrr::partial(label_both, sep = " = "))


#plt + 
#  scale_y_continuous(sec.axis = sec_axis(~ . , name = "SECOND Y AXIS", breaks = NULL, labels = NULL)) +
#  scale_x_continuous(sec.axis = sec_axis(~ . , name = "SECOND X AXIS", breaks = NULL, labels = NULL))
#+theme(axis.line.x.top = element_blank())
#+theme(axis.line.y.right = element_blank())