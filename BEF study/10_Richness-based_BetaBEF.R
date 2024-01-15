source('0_function.R')
###############################
folder = paste0('10_Richness')
dir.create(folder)
###############################
set.seed(123)
##################################################

# initialize project data
########################
load('project_16S_ASV')
Project_raw=get(paste0('project_',prefix))
# Env load
ENV <- sample_data(Project_raw)
Group = ENV$Elevation %>% as.factor(.)

# Examine effect of sample depth
Sample_Depth = 40000
SD_list = seq(5000,Sample_Depth,by=5000)
Depth = R2.a = P.a = P.b = R2.b = c()
for (i in 1:length(SD_list)) {
  Sample_Depth = SD_list[i]
  Project = Project_raw %>% rarefy_even_depth(.,Sample_Depth,rngseed=T)
  OTU <- otu_table(Project)
  OTU.dist = t(OTU) %>% vegdist(.)
  ENV <- sample_data(Project)
  
  ENV$Richness = specnumber(t(OTU))
  ENV$Dispersal = betadisper(OTU.dist,Elevation) %>% .$distance
  result.a = lm(Rs~Richness,data.frame(ENV)) %>% summary()
  result.b = lm(Rs~Dispersal,data.frame(ENV)) %>% summary()
  
  R2.a = c(R2.a,result.a$r.squared)
  P.a = c(P.a,result.a$coefficients['Richness','Pr(>|t|)'])
  
  R2.b = c(R2.b,result.b$r.squared)
  P.b = c(P.b,result.b$coefficients['Dispersal','Pr(>|t|)'])
  
  Depth = c(Depth,Sample_Depth)
}
Dt1 = data.frame(Depth,R2.a,R2.b,Data = 'Bacterial community')

########################
load('project_ITS_ASV')
Project_raw=get(paste0('project_',prefix))

Depth = R2.a = P.a = P.b = R2.b = c()
for (i in 1:length(SD_list)) {
  Sample_Depth = SD_list[i]
  Project = Project_raw %>% rarefy_even_depth(.,Sample_Depth,rngseed=T)
  OTU <- otu_table(Project)
  OTU.dist = t(OTU) %>% vegdist(.)
  ENV <- sample_data(Project)
  
  ENV$Richness = specnumber(t(OTU))
  ENV$Dispersal = betadisper(OTU.dist,Elevation) %>% .$distance
  result.a = lm(Rs~Richness,data.frame(ENV)) %>% summary()
  result.b = lm(Rs~Dispersal,data.frame(ENV)) %>% summary()
  
  R2.a = c(R2.a,result.a$r.squared)
  P.a = c(P.a,result.a$coefficients['Richness','Pr(>|t|)'])
  
  R2.b = c(R2.b,result.b$r.squared)
  P.b = c(P.b,result.b$coefficients['Dispersal','Pr(>|t|)'])
  
  Depth = c(Depth,Sample_Depth)
}
Dt2 = data.frame(Depth,R2.a,R2.b,Data = 'Fungal community')
Dt = rbind(Dt1,Dt2)

ggplot(Dt,aes(Depth,R2.a,color=Data))+geom_line()+geom_point(shape=21,size=3,fill='white')+theme_bw()+labs(x = "Depth of resampling reads", y = "R-squared of α-diversity",color ='Sequencing depth')+ylim(0,0.09)
ggsave(filename = paste0(folder,'/Sample_depth_on_alpha_bef.pdf'),dpi=600,width=12,height=10,units='cm')

ggplot(Dt,aes(Depth,R2.b,color=Data))+geom_line()+geom_point(shape=21,size=3,fill='white')+theme_bw()+labs(x = "Depth of resampling reads", y = "R-squared of β-diversity",color ='Sequencing depth')+ylim(0,0.09)
ggsave(filename = paste0(folder,'/Sample_depth_on_beta_bef.pdf'),dpi=600,width=12,height=10,units='cm')
