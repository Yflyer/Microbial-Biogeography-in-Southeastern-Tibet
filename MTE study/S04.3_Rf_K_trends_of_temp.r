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
dt = read.csv('S01/Phylo_Kgroup_RF-Richness_16S_ASV.csv')
Env = read.csv('S02/MTE_ENV.csv',row.names = 1)
Env = Env %>% select(Sample_ID,Rs,pH,TN,NH4,NO3,SOC,Temp.,Water_content,CN_ratio,HeE,Elevation)

Elevation = Env$Elevation

### Evaluate trends of Phylogeny and richness
plot_dt = dt  %>%  left_join(Env) %>% filter(K>3)
plot_dt$K = factor(plot_dt$K,levels = c(4:12))

##### MPD and Temp
stat_dt = plot_dt %>% group_by(K) %>% 
  summarise(Trend = formatC(lm_par(mpd,Temp.)[1], format = "e", digits = 3),
            R2 =format(lm_par(mpd,Temp.)[4],digits =3,nsmall=3),
            P = labelpstar(lm_par(mpd,Temp.)[2])) %>% 
  mutate(Text=paste0('K=',K,' E(a) = ',Trend,' R2 = ',R2,P))
stat_dt$Sign = 'Significant'
stat_dt$Sign[stat_dt$P==''] = 'Insignificant'
plot_dt2 = plot_dt %>% left_join(stat_dt)

ggplot(plot_dt2)+
  geom_point(alpha=0.5,aes(Temp.,mpd,color=K))+
  theme_bw()+
  guides(fill = 'none',linetype='none')+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        legend.position = 'right',
        legend.key.size = unit(0.5, 'cm'),
        legend.background = element_rect(colour = "black"),
        legend.text = element_text(size=12))+
  geom_smooth(method='lm',aes(Temp.,mpd,color=K,linetype=Sign),se=F)+
  scale_linetype_manual(values = c('dashed','solid'))+
  scale_color_discrete(name = "Slopes", labels = stat_dt$Text) +
  force_panelsizes(rows = unit(6, "cm"),cols = unit(6, "cm"))
ggsave(filename = paste0(folder,'/SI_Rf_Temp_mpd.pdf'),dpi=600,width=18,height=14,units='cm')

##### mntd and Temp
stat_dt = plot_dt %>% group_by(K) %>% 
  summarise(Trend = formatC(lm_par(mntd,Temp.)[1], format = "e", digits = 3),
            R2 =format(lm_par(mntd,Temp.)[4],digits =3,nsmall=3),
            P = labelpstar(lm_par(mntd,Temp.)[2])) %>% 
  mutate(Text=paste0('K=',K,' E(a) = ',Trend,' R2 = ',R2,P))
stat_dt$Sign = 'Significant'
stat_dt$Sign[stat_dt$P==''] = 'Insignificant'
plot_dt3 = plot_dt %>% left_join(stat_dt)

ggplot(plot_dt3)+
  geom_point(alpha=0.5,aes(Temp.,mntd,color=K))+
  theme_bw()+
  guides(fill = 'none',linetype='none')+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        legend.position = 'right',
        legend.key.size = unit(0.5, 'cm'),
        legend.background = element_rect(colour = "black"),
        legend.text = element_text(size=12))+
  geom_smooth(method='lm',aes(Temp.,mntd,color=K,linetype=Sign),se=F)+
  scale_linetype_manual(values = c('dashed','solid'))+
  scale_color_discrete(name = "Slopes", labels = stat_dt$Text) +
  force_panelsizes(rows = unit(6, "cm"),cols = unit(6, "cm"))
ggsave(filename = paste0(folder,'/SI_Rf_Temp_mntd.pdf'),dpi=600,width=18,height=14,units='cm')
