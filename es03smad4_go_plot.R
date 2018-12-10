library(ggplot2)
library(dplyr)
library(tidyr)

file = "/Users/jiejiaxu/kclab/SMAD4-RNA-seq/K.Chien_17_01-P8104/my_results/GOrilla/es03smad4_go_plot.txt"
data = read.delim(file,header = FALSE)
names(data) = c("GO_terms","FDR_qvalue","group")

colors = c("#B2B2FF","#FFB2B2","#B3F3B2")#brewer.pal(n=3, name="RdYlBu")
groups = c("day3_up","day6_up","day6_down")
for(k in seq(1,3)){
  data1 = filter(data,group == groups[k])
  data1 = data1[order(data1$FDR_qvalue),]
  data1$GO_terms = factor(data1$GO_terms,levels = rev(data1$GO_terms), ordered = T)
  p = ggplot(data = data1)+
    geom_bar(aes(x = GO_terms,y = -log10(FDR_qvalue)),fill = colors[k],stat = "identity",width = 0.5)+
    xlab("") + 
    theme_minimal()+ coord_flip()
  print(p)
}
