source('0_function.R')

#################### save dir #######################
folder = paste0('5_Robustness_analysis_TaxanomicResolution/')
dir.create(folder)
########################
model_result = data.frame(matrix(ncol = 5,nrow=7))
colnames(model_result) =c('Level','Mean.disper','Estimate','P_value','R2')

############### Env data pre
ENV <- read.csv('1_Data_Merge/Merge_dt.csv',row.names = 1,stringsAsFactors = T)
Elevation = ENV$Elevation %>% as.factor()

######################## Bac part
########################
load('project_16S_OTU')
project=get(paste0('project_',prefix))
OTU <- otu_table(project)
TAX <- tax_table(project) %>% data.frame()
OTU.dist = t(OTU) %>% vegdist()
ENV$Disper = betadisper(OTU.dist,Elevation) %>% .$distance
model_result$Mean.disper[1] = mean(ENV$Disper) %>% specify_decimal(.,3)
model_result$Se.disper[1] = sd(ENV$Disper,na.rm=T) %>% specify_decimal(.,3)

#### raw effect
result = lm(Rs~Disper,data=ENV) %>% summary() 
result.cof = result$coefficients %>% .['Disper',c('Estimate','Pr(>|t|)')]
model_result$Estimate[1]=result.cof['Estimate'] %>% specify_decimal(.,3)
model_result$P_value[1]=result.cof['Pr(>|t|)'] %>% specify_decimal(.,3)
model_result$R2[1]=result$adj.r.squared %>% specify_decimal(.,3)

model_result$Level[1]='OTU'
########################
load('project_16S_ASV')
project=get(paste0('project_',prefix))
OTU <- otu_table(project)
TAX <- tax_table(project) %>% data.frame()
OTU.dist = t(OTU) %>% vegdist()
ENV$Disper = betadisper(OTU.dist,Elevation) %>% .$distance
model_result$Mean.disper[2] = mean(ENV$Disper) %>% specify_decimal(.,3)
model_result$Se.disper[2] = sd(ENV$Disper,na.rm=T) %>% specify_decimal(.,3)

#### raw effect
result = lm(Rs~Disper,data=ENV) %>% summary() 
result.cof = result$coefficients %>% .['Disper',c('Estimate','Pr(>|t|)')]
model_result$Estimate[2]=result.cof['Estimate'] %>% specify_decimal(.,3)
model_result$P_value[2]=result.cof['Pr(>|t|)'] %>% specify_decimal(.,3)
model_result$R2[2]=result$adj.r.squared %>% specify_decimal(.,3)

model_result$Level[2]='ASV'

### Resolutions
Level = c("genus","family","order","class","phylum")
for (i in 1:length(Level)) {
  Taxa_dt <- TAX %>% unite("Taxa_index", kingdom:Level[i])
  OTU.dist = group_sum(t(OTU),Taxa_dt$Taxa_index) %>% vegdist() 
  
  ENV$Disper = betadisper(OTU.dist,Elevation) %>% .$distance
  model_result$Mean.disper[i+2] = mean(ENV$Disper) %>% specify_decimal(.,3)
  model_result$Se.disper[i+2] = sd(ENV$Disper,na.rm=T) %>% specify_decimal(.,3)
  #### raw effect
  result = lm(Rs~Disper,data=ENV) %>% summary() 
  result.cof = result$coefficients %>% .['Disper',c('Estimate','Pr(>|t|)')]
  model_result$Estimate[i+2]=result.cof['Estimate'] %>% specify_decimal(.,3)
  model_result$P_value[i+2]=result.cof['Pr(>|t|)'] %>% specify_decimal(.,3)
  model_result$R2[i+2]=result$adj.r.squared %>% specify_decimal(.,3)
  model_result$Level[i+2]=Level[i]
}

write.csv(model_result,paste0(folder,'/Bac_level_effects_on_BEF.csv'))

######################## Fug part
########################
load('project_ITS_OTU')
project=get(paste0('project_',prefix))
OTU <- otu_table(project)
TAX <- tax_table(project) %>% data.frame()
OTU.dist = t(OTU) %>% vegdist()
ENV$Disper = betadisper(OTU.dist,Elevation) %>% .$distance
model_result$Mean.disper[1] = mean(ENV$Disper) %>% specify_decimal(.,3)
model_result$Se.disper[1] = sd(ENV$Disper,na.rm=T) %>% specify_decimal(.,3)

#### raw effect
result = lm(Rs~Disper,data=ENV) %>% summary() 
result.cof = result$coefficients %>% .['Disper',c('Estimate','Pr(>|t|)')]
model_result$Estimate[1]=result.cof['Estimate'] %>% specify_decimal(.,3)
model_result$P_value[1]=result.cof['Pr(>|t|)'] %>% specify_decimal(.,3)
model_result$R2[1]=result$adj.r.squared %>% specify_decimal(.,3)

model_result$Level[1]='OTU'
########################
load('project_ITS_ASV')
project=get(paste0('project_',prefix))
OTU <- otu_table(project)
TAX <- tax_table(project) %>% data.frame()
OTU.dist = t(OTU) %>% vegdist()
ENV$Disper = betadisper(OTU.dist,Elevation) %>% .$distance
model_result$Mean.disper[2] = mean(ENV$Disper) %>% specify_decimal(.,3)
model_result$Se.disper[2] = sd(ENV$Disper,na.rm=T) %>% specify_decimal(.,3)

#### raw effect
result = lm(Rs~Disper,data=ENV) %>% summary() 
result.cof = result$coefficients %>% .['Disper',c('Estimate','Pr(>|t|)')]
model_result$Estimate[2]=result.cof['Estimate'] %>% specify_decimal(.,3)
model_result$P_value[2]=result.cof['Pr(>|t|)'] %>% specify_decimal(.,3)
model_result$R2[2]=result$adj.r.squared %>% specify_decimal(.,3)

model_result$Level[2]='ASV'
### Resolutions
Level = c("genus","family","order","class","phylum")
for (i in 1:length(Level)) {
  Taxa_dt <- TAX %>% unite("Taxa_index", kingdom:Level[i])
  OTU.dist = group_sum(t(OTU),Taxa_dt$Taxa_index) %>% vegdist() 
  
  ENV$Disper = betadisper(OTU.dist,Elevation) %>% .$distance
  model_result$Mean.disper[i+2] = mean(ENV$Disper) %>% specify_decimal(.,3)
  model_result$Se.disper[i+2] = sd(ENV$Disper,na.rm=T) %>% specify_decimal(.,3)
  #### raw effect
  result = lm(Rs~Disper,data=ENV) %>% summary() 
  result.cof = result$coefficients %>% .['Disper',c('Estimate','Pr(>|t|)')]
  model_result$Estimate[i+2]=result.cof['Estimate'] %>% specify_decimal(.,3)
  model_result$P_value[i+2]=result.cof['Pr(>|t|)'] %>% specify_decimal(.,3)
  model_result$R2[i+2]=result$adj.r.squared %>% specify_decimal(.,3)
  model_result$Level[i+2]=Level[i]
}

write.csv(model_result,paste0(folder,'/Fug_level_effects_on_BEF.csv'))