source('0_function.R')
library(ieggr)
library(ggh4x)
#################### save dir #######################
folder = paste0('S04')
dir.create(folder)
########### input ##############
# set input
load('project_16S_ASV')
project = get(paste0('project_',prefix))
###
ENV = read.csv('S02/MTE_ENV.csv',row.names = 1,stringsAsFactors = T) 
Elevation = ENV$Elevation

dt = read.csv('S01/Phylo_Kgroup_LocalScale_16S_ASV_layers.csv',row.names = 1) %>% mutate(Rf='No-detrending') %>% left_join(ENV)
dt$K = factor(dt$K,levels = c(3:12))

### Trends of richness by K
stat_dt = dt %>% group_by(K) %>% 
  summarise(Trend = format(lm_par(log10(Richness),Inv.kT)[1],digits=3,nsmall=3),
            R2 =format(lm_par(log10(Richness),Inv.kT)[4],digits =3,nsmall=3),
            P = labelpstar(lm_par(log10(Richness),Inv.kT)[2])) %>% 
  mutate(Text=paste0('K=',K,' E(a) = ',Trend,' R2 = ',R2,P))

ggplot(dt)+
  geom_point(alpha=0.5,aes(Temp.,log10(Richness),color=K))+
  theme_bw()+
  guides(fill = 'none')+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        legend.position = 'right',
        legend.key.size = unit(0.5, 'cm'),
        legend.background = element_rect(colour = "black"),
        legend.text = element_text(size=12))+
  geom_smooth(method='lm',aes(Temp.,log10(Richness),color=K))+
  scale_color_discrete(name = "Activate energys of MTE", labels = stat_dt$Text) +
  #labs(x = "Occurence of elevations (K)", y = "Phylogenetic distance of nearest taxon",color ='Detrending')+
  #facet_wrap(~Rf)+
  #ylim(4.5, 8.5)+
  force_panelsizes(rows = unit(6, "cm"),cols = unit(6, "cm"))
ggsave(filename = paste0(folder,'/Temp_K_MTE.pdf'),dpi=600,width=16,height=14,units='cm')

##### Richness proportions change
dt = dt %>% 
  group_by(Sample_ID) %>% 
  mutate(Total_richness= sum(Richness),
         Total_reads = sum(Reads.sum),
         Richness_prop = Richness/Total_richness*100,
         Reads_prop = Reads.sum/Richness_prop*100)
### stat
stat_dt = dt %>% group_by(K) %>% 
  summarise(Trend = format(lm_par(Richness_prop,Temp.)[1],digits=3,nsmall=3),
            R2 =format(lm_par(Richness_prop,Temp.)[4],digits =3,nsmall=3),
            P = labelpstar(lm_par(Richness_prop,Temp.)[2])) %>% 
  mutate(Text=paste0('K=',K,' Slope = ',Trend,'% R2 = ',R2,P))

ggplot(dt,aes(x=Temp.,y=Richness_prop)) + 
  geom_point(aes(color=K),alpha=0.5)+
  geom_smooth(aes(color=K),method = 'lm')+
  theme_bw()+
  guides(fill = 'none')+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        legend.position = 'right',
        legend.key.size = unit(0.5, 'cm'),
        legend.background = element_rect(colour = "black"),
        legend.text = element_text(size=12))+
  scale_color_discrete(name = "Proportional change of richness", labels = stat_dt$Text) +
  scale_y_continuous(name = "Proportion of richness",n.breaks = 5,breaks = c(0,10,20,30),labels = c('0%','10%','20%','30%'))+
  #ylab('Proportion of richness')+
  #ylim(4.5, 8.5)+
  force_panelsizes(rows = unit(6, "cm"),cols = unit(6, "cm"))
ggsave(filename = paste0(folder,'/Temp_K_Richness_proportion.pdf'),dpi=600,width=16,height=14,units='cm')