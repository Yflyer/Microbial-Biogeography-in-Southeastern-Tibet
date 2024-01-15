source('0_function.R')
library(ggh4x)
#################### save dir #######################
folder = paste0('S01')
dir.create(folder)
######################
### regression plot
lm_par = function(Richness,Temp.){
  decay <- lm(Richness~Temp.)
  lm.slope = decay$coefficients ['Temp.'] %>% round(.,3)
  lm.intercept = decay$coefficients ['(Intercept)'] %>% round(.,3)
  lm.P = lmp(decay) %>% labelpstar()
  R2= summary(decay)$adj.r.squared %>% round(.,3)
  N = length(summary(decay)$residuals)
  return(c(lm.slope,lm.P,lm.intercept,R2,N)
  )
}
####################
### theme set
PTs = theme(axis.text.y   = element_text(size=12),
            axis.text.x   = element_text(size=12),
            axis.title.y  = element_text(size=14),
            axis.title.x  = element_text(size=14),
            panel.background=element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1.5))
####################
Env = read.csv('S01/MTE_ENV.csv',row.names = 1)

####
load('project_16S_ASV')
project = get(paste0('project_',prefix))
OTU = otu_table(project) %>% data.frame()
plot_data = sample_data(project) %>% data.frame %>% 
  mutate(Sample_ID = sample_names(project)) %>% 
  select ( Sample_ID,Richness,Shannon,Invsimpson, PD) %>% full_join(Env)

# bacteria richness
Fit_par = lm_par(plot_data$Richness,plot_data$Temp.)
lm_txt=paste0('y = ',Fit_par[1],' x +',Fit_par[3],'\n','R2 = ',Fit_par[4],Fit_par[2],' N = ',Fit_par[5])

ggplot(plot_data,aes(x=Temp.,y=Richness))+geom_point()+
  geom_smooth(method='lm',color='black',se=F,linewidth=1.5)+
  annotate("text",x = 20,y=11000,label = lm_txt,size=3,vjust=1,hjust=0)+
  labs(x = "Soil temperature (°C)", y = "Richness")+
  PTs+
  force_panelsizes(rows = unit(5, "cm"),cols = unit(7, "cm"))

ggsave(filename = paste0(folder,'/Temperature_dependency_Soil_Bac_Richness.pdf'),dpi=900,width=12,height=8,units='cm')

# bacteria shannon
Fit_par = lm_par(plot_data$Shannon,plot_data$Temp.)
lm_txt=paste0('y = ',Fit_par[1],' x +',Fit_par[3],'\n','R2 = ',Fit_par[4],Fit_par[2],' N = ',Fit_par[5])

ggplot(plot_data,aes(x=Temp.,y=Shannon))+geom_point()+
    geom_smooth(method='lm',color='black',se=F,size=1.5)+
    annotate("text",x = 20,y=9,label = lm_txt,size=3,vjust=1,hjust=0)+
  labs(x = "Soil temperature (°C)", y = "Shannon index")+
  PTs+
  force_panelsizes(rows = unit(5, "cm"),cols = unit(7, "cm"))

ggsave(filename = paste0(folder,'/Temperature_dependency_Soil_Bac_shannon.pdf'),dpi=900,width=12,height=8,units='cm')

# bacteria invsimpson
Fit_par = lm_par(plot_data$Invsimpson,plot_data$Temp.)
lm_txt=paste0('y = ',Fit_par[1],' x +',Fit_par[3],'\n','R2 = ',Fit_par[4],Fit_par[2],' N = ',Fit_par[5])

ggplot(plot_data,aes(x=Temp.,y=Invsimpson))+geom_point()+
  geom_smooth(method='lm',color='black',se=F,size=1.5)+
  annotate("text",x = 20,y=2400,label = lm_txt,size=3,vjust=1,hjust=0)+
  labs(x = "Soil temperature (°C)", y = "1/(Simpson index)")+
  PTs+
  force_panelsizes(rows = unit(5, "cm"),cols = unit(7, "cm"))

ggsave(filename = paste0(folder,'/Temperature_dependency_Soil_Bac_invsimpson.pdf'),dpi=900,width=12,height=8,units='cm')

# bacteria PD
Fit_par = lm_par(plot_data$PD,plot_data$Temp.)
lm_txt=paste0('y = ',Fit_par[1],' x +',Fit_par[3],'\n','R2 = ',Fit_par[4],Fit_par[2],' N = ',Fit_par[5])
#lm(PD~Temp.,plot_data) %>% summary()
ggplot(plot_data,aes(x=Temp.,y=PD))+geom_point()+
  geom_smooth(method='lm',color='black',se=F,size=1.5)+
  annotate("text",x = 20,y=330,label = lm_txt,size=3,vjust=1,hjust=0)+
  labs(x = "Soil temperature (°C)", y = "Faith's PD")+
  PTs+
  force_panelsizes(rows = unit(5, "cm"),cols = unit(7, "cm"))
ggsave(filename = paste0(folder,'/Temperature_dependency_Soil_Bac_PD.pdf'),dpi=900,width=12,height=8,units='cm')

########################################################################
# fungal community
load('project_ITS_ASV')
project = get(paste0('project_',prefix))
OTU = otu_table(project) %>% data.frame()
plot_data = sample_data(project) %>% data.frame %>% 
  mutate(Sample_ID = sample_names(project)) %>% 
  select ( Sample_ID,Richness,Shannon,Invsimpson) %>% full_join(Env)

#  richness
Fit_par = lm_par(plot_data$Richness,plot_data$Temp.)
lm_txt=paste0('y = ',Fit_par[1],' x +',Fit_par[3],'\n','R2 = ',Fit_par[4],Fit_par[2],' N = ',Fit_par[5])

ggplot(plot_data,aes(x=Temp.,y=Richness))+geom_point()+
  geom_smooth(method='lm',color='black',se=F,size=1.5)+
  annotate("text",x = 20,y=1800,label = lm_txt,size=3,vjust=1,hjust=0)+
  labs(x = "Soil temperature (°C)", y = "Richness")+
  PTs+
  force_panelsizes(rows = unit(5, "cm"),cols = unit(7, "cm"))

ggsave(filename = paste0(folder,'/Temperature_dependency_Soil_Fug_Richness.pdf'),dpi=900,width=12,height=8,units='cm')

# bacteria shannon
Fit_par = lm_par(plot_data$Shannon,plot_data$Temp.)
lm_txt=paste0('y = ',Fit_par[1],' x +',Fit_par[3],'\n','R2 = ',Fit_par[4],Fit_par[2],' N = ',Fit_par[5])

ggplot(plot_data,aes(x=Temp.,y=Shannon))+geom_point()+
  geom_smooth(method='lm',color='black',se=F,size=1.5)+
  annotate("text",x = 20,y=6.5,label = lm_txt,size=3,vjust=1,hjust=0)+
  labs(x = "Soil temperature (°C)", y = "Shannon index")+
  PTs+
  force_panelsizes(rows = unit(5, "cm"),cols = unit(7, "cm"))

ggsave(filename = paste0(folder,'/Temperature_dependency_Soil_Fug_shannon.pdf'),dpi=900,width=12,height=8,units='cm')

# bacteria Fug_invsimpson
Fit_par = lm_par(plot_data$Invsimpson,plot_data$Temp.)
lm_txt=paste0('y = ',Fit_par[1],' x +',Fit_par[3],'\n','R2 = ',Fit_par[4],Fit_par[2],' N = ',Fit_par[5])

ggplot(plot_data,aes(x=Temp.,y=Invsimpson))+geom_point()+
  geom_smooth(method='lm',color='black',se=F,size=1.5)+
  annotate("text",x = 20,y=105,label = lm_txt,size=3,vjust=1,hjust=0)+
  labs(x = "Soil temperature (°C)", y = "1/(Simpson index)")+
  PTs+
  force_panelsizes(rows = unit(5, "cm"),cols = unit(7, "cm"))

ggsave(filename = paste0(folder,'/Temperature_dependency_Soil_Fug_invsimpson.pdf'),dpi=900,width=12,height=8,units='cm')

