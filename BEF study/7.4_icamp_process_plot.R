source('0_function.R')
library(ieggr)
library(tidyr)

#################### save dir #######################
folder = paste0('7_icamp')
dir.create(folder)
note ='all' # for optional subset
##################################################

####################
PTs <- theme(text = element_text(size=16),
             axis.text.x = element_text(angle = 45,hjust=1),
             #axis.ticks = element_blank(),
             panel.background=element_blank(),
             #legend.position = "none",
             panel.border = element_rect(colour = "black", fill=NA, size=1))

plot_data = read.csv('7_icamp_16S_ASV/16S_ASV.ProcessImportance_EachGroup.csv')
load('project_16S_ASV')
project = get(paste0('project_',prefix))

d1 = plot_data[1:12,] %>% 
  group_by(Method) %>% 
  summarise(`Heterogeneous selection` = mean(HeS),
            `Homogeneous selection` = mean(HoS),
            `Dispersal limitation` = mean(DL),
            `Drift` = mean(DR),
            `Homogenizing dispersal` = mean(HD)) %>% 
  mutate(Dataset = prefix, Group = 'in each elevation')
d2 = plot_data[13:78,] %>% 
  group_by(Method) %>% 
  summarise(`Heterogeneous selection` = mean(HeS),
            `Homogeneous selection` = mean(HoS),
            `Dispersal limitation` = mean(DL),
            `Drift` = mean(DR),
            `Homogenizing dispersal` = mean(HD)) %>% 
  mutate(Dataset = prefix, Group = 'between elevations')

plot_data = read.csv('7_icamp_ITS_ASV/ITS_ASV.ProcessImportance_EachGroup.csv')
load('project_ITS_ASV')
project = get(paste0('project_',prefix))

d3 = plot_data[1:12,] %>% 
  group_by(Method) %>% 
  summarise(`Heterogeneous selection` = mean(HeS),
            `Homogeneous selection` = mean(HoS),
            `Dispersal limitation` = mean(DL),
            `Drift` = mean(DR),
            `Homogenizing dispersal` = mean(HD)) %>% 
  mutate(Dataset = prefix, Group = 'in each elevation')
d4 = plot_data[13:78,] %>% 
  group_by(Method) %>% 
  summarise(`Heterogeneous selection` = mean(HeS),
            `Homogeneous selection` = mean(HoS),
            `Dispersal limitation` = mean(DL),
            `Drift` = mean(DR),
            `Homogenizing dispersal` = mean(HD)) %>% 
  mutate(Dataset = prefix, Group = 'between elevations')

bar_dt = Reduce(rbind,list(d1,d2,d3,d4)) %>% 
  gather(key = 'Process',value='Ratio',-Method,-Group,-Dataset) %>%
  mutate(Group = factor(Group,levels=c('in each elevation','between elevations')),
         Process = factor(Process,levels=c('Heterogeneous selection',
                                          'Homogeneous selection',
                                          'Homogenizing dispersal', 
                                           'Dispersal limitation',
                                          'Drift')))


library(ggplot2)
library(ggh4x)
ggplot(bar_dt) + 
  geom_bar(aes(Process,Ratio,fill=Group),stat = 'identity',position = 'dodge')+
  PTs +
  labs(x = "Stochastic processes", y = "Percentage")+
  geom_text(x = 1.5, y = 12, label = "***")+
  force_panelsizes(rows = unit(6, "cm"),cols = unit(8, "cm"))+
  facet_grid(~Dataset)

