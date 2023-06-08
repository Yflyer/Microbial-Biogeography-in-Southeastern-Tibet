# AO-plot
library(ggh4x)
source('0_function.R')
#################### save dir #######################
folder = paste0('8_A&O_analysis')
dir.create(folder)

##############################
lm_par = function(Y,X){
  decay <- lm(Y~X)
  lm.slope = decay$coefficients ['X'] %>% round(.,3)
  lm.intercept = decay$coefficients ['(Intercept)'] %>% round(.,3)
  lm.P = lmp(decay) %>% labelpstar()
  R2= summary(decay)$adj.r.squared %>% round(.,3)
  N = length(summary(decay)$residuals)
  return(c(lm.slope,lm.P,lm.intercept,R2,N))
}

############### Env data pre
ENV <- read.csv('1_Data_Merge/Merge_dt.csv',row.names = 1,stringsAsFactors = T)
Elevation = ENV$Elevation %>% as.factor()
################### Fug ############################
### load project  
load('project_ITS_ASV')
project = get(paste0('project_',prefix)) 
OTU <- otu_table(project)

load("7_stochasticity/ITS_ASV_all_st_boot")
Dt_st_group = get("ITS_ASV_all_st_boot") %>% .$summary %>% filter(Index == 'NST') %>% select(Plot = Group,NST = obs,NST_sd = stdev)

######################################################
#### create group occur data
count_table=group_count(table = OTU, group = Elevation) %>% as.data.frame(.)
group_occur= rowSums(count_table>0)
Dt = data.frame()
for (i in 1:ncol(count_table)) {
  Species_proportion = apply(count_table,2,function(x)sum(group_occur[x>0]==i)/sum(x>0))
  Sub_dt = data.frame(Elevation = as.integer(levels(Elevation)),Species_proportion,Species_occupancy=i)
  Dt =rbind(Dt,Sub_dt)
}
Dt =  data.frame(ENV[,c('Plot','Elevation')]) %>% distinct(.) %>% left_join(Dt,.) %>% left_join(.,Dt_st_group) 
Dt$Prefix = 'Fungal community'

AO_dt_group = Dt

################### Bac ############################
### load project  
load('project_16S_ASV')
project = get(paste0('project_',prefix)) 
OTU <- otu_table(project)
ENV <- sample_data(project)

load("7_stochasticity/16S_ASV_all_st_boot")
Dt_st_group = get("16S_ASV_all_st_boot") %>% .$summary %>% filter(Index == 'NST') %>% select(Plot = Group,NST = obs,NST_sd = stdev)

#############################################################
#### create group occur data
count_table=group_count(table = OTU, group = Elevation) %>% as.data.frame(.)
group_occur= rowSums(count_table>0)
Dt = data.frame()
for (i in 1:ncol(count_table)) {
  Species_proportion = apply(count_table,2,function(x)sum(group_occur[x>0]==i)/sum(x>0))
  Sub_dt = data.frame(Elevation = as.integer(levels(Elevation)),Species_proportion,Species_occupancy=i)
  Dt =rbind(Dt,Sub_dt)
}
Dt =  data.frame(ENV[,c('Plot','Elevation')]) %>% distinct(.) %>% left_join(Dt,.) %>% left_join(.,Dt_st_group)
Dt$Prefix = 'Bacterial community'

AO_dt_group = rbind(AO_dt_group,Dt)
AO_dt_group$Elevation = as.factor(AO_dt_group$Elevation )
AO_dt_group$Species_occupancy = factor(AO_dt_group$Species_occupancy,levels = 12:1)

#################### Color scheme
slm = scale_linetype_manual(values = c('solid','dashed'))
sfm = scale_fill_manual(values = colorRampPalette(colors = c('#ec7696','#f1c4cd','#d8e3e7','#d0dfe6','#11659a'))(12))
scm = scale_color_manual(values = colorRampPalette(colors = c('#ec7696','#f1c4cd','#d8e3e7','#d0dfe6','#11659a'))(12))
ssm = scale_shape_manual(values = c(16,21))

####################
# AO pattern
ggplot(AO_dt_group)+
  geom_bar(aes(Elevation,Species_proportion,fill=Species_occupancy),stat = 'identity',width=.9)+
  sfm+
  labs(y = "Species proportion(%)", x = "Elevation(m)")+
  scale_y_continuous(labels = scales::percent,limits = c(0,1.15),breaks = seq(0,1,0.5))+
  facet_grid(Prefix~.) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5),
        panel.grid = element_blank(),
        panel.spacing = unit(0.2, "lines"),
        legend.key.size = unit(0.1, "cm"),
        legend.key.width = unit(0.5, "cm"))+
  force_panelsizes(rows = unit(2.5, "cm"),cols = unit(6.5, "cm"))
ggsave(filename = paste0(folder,'/Elevation_AO.pdf'),dpi=300,width=14,height=9,units='cm')

##################  NST trends
plot_dt = AO_dt_group %>%
  distinct(Elevation,NST,NST_sd,Prefix)
ggplot(plot_dt,aes(Elevation,NST))+geom_col(fill='grey80', color="black",width=.9)+
  scale_color_viridis()+
  geom_errorbar(aes(ymin=NST, ymax=NST+NST_sd), width=.2,
                position=position_dodge(.9))+
  scale_y_continuous(labels = scales::percent,limits = c(0,1.15),breaks = seq(0,1,0.5))+
  labs(y = "Stochasticity(%)", x = "Elevation(m)")+
  facet_grid(Prefix~.) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5),
        panel.grid = element_blank(),
        panel.spacing = unit(0.2, "lines"))+
  force_panelsizes(rows = unit(2.5, "cm"),cols = unit(6.5, "cm"))
ggsave(filename = paste0(folder,'/Elevation_NST.pdf'),dpi=300,width=14,height=9,units='cm')

######################################################
################## Trend stat data
lm_dt = AO_dt_group %>%
  group_by(Prefix,Species_occupancy) %>%
  summarise(Slope=lm_par(Species_proportion,NST)[1],
            P = lm_par(Species_proportion,NST)[2],
            R2=lm_par(Species_proportion,NST)[4],
            Text = paste0('y = ',lm_par(Species_proportion,NST)[1],'x + ',lm_par(Species_proportion,NST)[3],'   R2 = ',lm_par(Species_proportion,NST)[4],lm_par(Species_proportion,NST)[2]))
lm_dt$P [lm_dt$P %in% c('')] = 'Unsignif.'
lm_dt$P [lm_dt$P != 'Unsignif.'] = 'Signif.'

##################  bac
plot_dt = AO_dt_group %>%
  filter(Prefix=='Bacterial community') %>%
  left_join(lm_dt)
ggplot(plot_dt)+
  geom_smooth(aes(NST,Species_proportion,group=Species_occupancy, color=Species_occupancy,linetype=P),method = 'lm',se=F)+scm+
  labs(y = "Species proportion(%)", x = "Stochasticity(%)")+
  scale_x_continuous(labels = scales::percent,limits =c(0.6,NA))+
  scale_y_continuous(labels = scales::percent,limits =c(0,0.17))+
  facet_grid(~Prefix) +
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'null')+
  force_panelsizes(rows = unit(4, "cm"),cols = unit(3, "cm"))
ggsave(filename = paste0(folder,'/NST_AO_bacterial_series.pdf'),dpi=300,width=5,height=6,units='cm')

##################  fug
plot_dt = AO_dt_group %>%
  filter(Prefix=='Fungal community') %>%
  left_join(lm_dt)
ggplot(plot_dt)+
  geom_smooth(aes(NST,Species_proportion,group=Species_occupancy, color=Species_occupancy,linetype=P),method = 'lm',se=F)+scm+
  labs(y = "Species proportion(%)", x = "Stochasticity(%)")+
  labs(y = "Species proportion(%)", x = "Stochasticity(%)")+
  scale_y_continuous(labels = scales::percent,limits =c(0,NA))+
  scale_x_continuous(labels = scales::percent,limits =c(0.1,1))+
  facet_grid(~Prefix) +
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'null')+
  force_panelsizes(rows = unit(4, "cm"),cols = unit(3, "cm"))
ggsave(filename = paste0(folder,'/NST_AO_fungal_series.pdf.pdf'),dpi=300,width=5,height=6,units='cm')
  
