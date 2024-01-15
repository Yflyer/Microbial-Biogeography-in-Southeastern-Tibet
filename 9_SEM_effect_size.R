source('0_function.R')
library(lavaan)
library(semTools)
library(ggh4x)
####
#################### save dir #######################
folder = paste0('9_sem')
dir.create(folder)
#################### load data #######################
plot_dt <- read.csv('9_sem/Effect size.csv',row.names = 1,stringsAsFactors = T)
colnames(plot_dt)
ggplot(plot_dt)+
  geom_col(aes(x=Diversity,y=Effect.size,fill=Type,linetype=Significance),color='black')+
  scale_linetype_manual(values = c('dashed','solid'))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45))
ggsave(filename = paste0(folder,'/Effect.size.pdf'),dpi=600,width=12,height=10,units='cm')
