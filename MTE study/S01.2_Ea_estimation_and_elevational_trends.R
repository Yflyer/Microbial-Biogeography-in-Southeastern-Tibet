source('0_function.R')
library(ggh4x)
library(rsq)
#################### save dir #######################
folder = paste0('S01')
dir.create(folder)
####################
### regression plot
lm_par = function(Richness,Inv.kT){
  decay <- lm(Richness~Inv.kT)
  summ = summary(decay)$coefficients
  lm.slope = summ[2,1] %>% round(.,3)
  lm.slope.sd = summ[2,2] %>% round(.,3)
  lm.intercept = summ[1,1] %>% round(.,3)
  lm.P = lmp(decay) #%>% labelpstar()
  R2= summary(decay)$adj.r.squared %>% round(.,3)
  return(c(lm.slope,lm.P,lm.intercept,R2,lm.slope.sd))
}

####################
plot_data = read.csv('S01/MTE_ENV.csv',row.names = 1) #%>% filter(Elevation!=2013)

Elevation=plot_data$Elevation
plot_data$Elevation = as.factor(plot_data$Elevation)

dt = data.frame()
#################### chao1 estimation
# 100% Identity
load('project_16S_ASV')
project = get(paste0('project_16S_ASV')) #%>% subset_samples(Elevation!=2013)
OTU <- otu_table(project) %>% data.frame() %>% .[,plot_data$Sample_ID]
plot_data$Soil.Bac.Chao1 = estimateR(t(OTU))['S.chao1',]
plot_data$Soil.Bac.Richness = specnumber(t(OTU))
# Ea
Fit_par = lm_par(log(plot_data$Soil.Bac.Chao1),plot_data$Inv.kT)
Ea = Fit_par[1] %>% round(.,3)
# Temp
Fit_par = lm_par(plot_data$Soil.Bac.Richness,plot_data$Temp.)
Eestimate.temp = Fit_par[1] %>% round(.,3)
Eestimate.sd.temp=Fit_par[5] %>% round(.,3)
P.temp = Fit_par[2] %>% round(.,3)
R2.temp = Fit_par[4] %>% round(.,3)
# Ele
Fit_par = lm_par(plot_data$Soil.Bac.Richness,as.numeric(plot_data$Elevation))
Eestimate.ele = Fit_par[1] %>% round(.,3)
Eestimate.sd.ele=Fit_par[5] %>% round(.,3)
P.ele = Fit_par[2] %>% round(.,3)
R2.ele = Fit_par[4] %>% round(.,3)
# multi
formula.control = paste0('Soil.Bac.Richness~HeE+Temp.+Water_content+pH+SOC+TN+CN_ratio+NH4') %>% as.formula()
mlm <- lm(formula.control, data = plot_data)
slm = step(mlm,direction = "backward")
sum.slm = summary(slm)
rsq.result = rsq.partial(slm,adj = T)
P.R2.prop = (rsq.result$partial.rsq[2] / sum.slm$adj.r.squared) %>% round(.,3)

Result = c(Eestimate.ele,Eestimate.sd.ele,P.ele,R2.ele,
           Eestimate.temp,Eestimate.sd.temp,P.temp,R2.temp,
           P.R2.prop,Ea)
dt=rbind(dt,Result)

# 99% Identity
load('Project_16S_99')
project = get(paste0('Project_16S_99')) %>% rarefy_even_depth(sample.size=40000)
OTU <- otu_table(project) %>% data.frame() %>% .[,plot_data$Sample_ID]
plot_data$Soil.Bac.Chao1 = estimateR(t(OTU))['S.chao1',]
plot_data$Soil.Bac.Richness = specnumber(t(OTU))
# Ea
Fit_par = lm_par(log(plot_data$Soil.Bac.Chao1),plot_data$Inv.kT)
Ea = Fit_par[1] %>% round(.,3)
# Temp
Fit_par = lm_par(plot_data$Soil.Bac.Richness,plot_data$Temp.)
Eestimate.temp = Fit_par[1] %>% round(.,3)
Eestimate.sd.temp=Fit_par[5] %>% round(.,3)
P.temp = Fit_par[2] %>% round(.,3)
R2.temp = Fit_par[4] %>% round(.,3)
# Ele
Fit_par = lm_par(plot_data$Soil.Bac.Richness,as.numeric(plot_data$Elevation))
Eestimate.ele = Fit_par[1] %>% round(.,3)
Eestimate.sd.ele=Fit_par[5] %>% round(.,3)
P.ele = Fit_par[2] %>% round(.,3)
R2.ele = Fit_par[4] %>% round(.,3)
# multi
formula.control = paste0('Soil.Bac.Richness~HeE+Temp.+Water_content+pH+SOC+TN+CN_ratio+NH4') %>% as.formula()
mlm <- lm(formula.control, data = plot_data)
slm = step(mlm,direction = "backward")
sum.slm = summary(slm)
rsq.result = rsq.partial(slm,adj = T)
P.R2.prop = (rsq.result$partial.rsq[2] / sum.slm$adj.r.squared) %>% round(.,3)

Result = c(Eestimate.ele,Eestimate.sd.ele,P.ele,R2.ele,
  Eestimate.temp,Eestimate.sd.temp,P.temp,R2.temp,
  P.R2.prop,Ea)
dt=rbind(dt,Result)

# 97% Identity
load('project_16S_OTU')
project = get(paste0('project_16S_OTU'))
OTU <- otu_table(project) %>% data.frame() %>% .[,plot_data$Sample_ID]
plot_data$Soil.Bac.Chao1 = estimateR(t(OTU))['S.chao1',]
plot_data$Soil.Bac.Richness = specnumber(t(OTU))
# Ea
Fit_par = lm_par(log(plot_data$Soil.Bac.Chao1),plot_data$Inv.kT)
Ea = Fit_par[1] %>% round(.,3)
# Temp
Fit_par = lm_par(plot_data$Soil.Bac.Richness,plot_data$Temp.)
Eestimate.temp = Fit_par[1] %>% round(.,3)
Eestimate.sd.temp=Fit_par[5] %>% round(.,3)
P.temp = Fit_par[2] %>% round(.,3)
R2.temp = Fit_par[4] %>% round(.,3)
# Ele
Fit_par = lm_par(plot_data$Soil.Bac.Richness,as.numeric(plot_data$Elevation))
Eestimate.ele = Fit_par[1] %>% round(.,3)
Eestimate.sd.ele=Fit_par[5] %>% round(.,3)
P.ele = Fit_par[2] %>% round(.,3)
R2.ele = Fit_par[4] %>% round(.,3)
# multi
formula.control = paste0('Soil.Bac.Richness~HeE+Temp.+Water_content+pH+SOC+TN+CN_ratio+NH4') %>% as.formula()
mlm <- lm(formula.control, data = plot_data)
slm = step(mlm,direction = "backward")
sum.slm = summary(slm)
rsq.result = rsq.partial(slm,adj = T)
P.R2.prop = (rsq.result$partial.rsq[2] / sum.slm$adj.r.squared) %>% round(.,3)

### Combine result
Result = c(Eestimate.ele,Eestimate.sd.ele,P.ele,R2.ele,
           Eestimate.temp,Eestimate.sd.temp,P.temp,R2.temp,
           P.R2.prop,Ea)
dt=rbind(dt,Result)
colnames(dt)= c('Eestimate.ele','Eestimate.sd.ele','P.ele','R2.ele',
                'Eestimate.temp','Eestimate.sd.temp','P.temp','R2.temp',
                'P.R2.prop','Ea')
dt$Data = c('ASV','OTU(99%)','OTU(97%)')
write.csv(dt,file = paste0(folder,'/MTE_stat_overall.csv'))

