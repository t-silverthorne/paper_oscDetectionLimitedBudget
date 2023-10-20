c(library(ggplotify),
library(ggplot2),
library(patchwork))

df=data.frame(x=1,y=1)

pa=ggplot(df,aes(x=x,y=y))+geom_point()
pb=pa
pc=pa

p1 = pa/pb

p1 = p1+plot_annotation(title = 'test1')
p1=as.ggplot(p1)
p2 = p1

(p1/p2) + plot_annotation(tag_levels='A')
