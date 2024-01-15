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
plot_data$Elevation = factor(plot_data$Elevation,sort(unique(plot_data$Elevation)))

dt = data.frame()
#################### Fitting type
Variable = Fitting_type = Data_type = AIC = R2 = P_value = c()
# 100% Identity
load('project_16S_ASV')
project = get(paste0('project_16S_ASV')) #%>% subset_samples(Elevation!=2013)
OTU <- otu_table(project) %>% data.frame() %>% .[,plot_data$Sample_ID]

### Temp.
# linear
plot_data$Richness = specnumber(t(OTU))
mod = lm(Richness~Temp.,plot_data) 
AIC = c(AIC,AIC(mod)) %>% round(3)
R2 = c(R2,rsq(mod,adj = T)) %>% round(3)
P_value = c(P_value,lmp(mod)) %>% round(3)
Variable = c(Variable,'Soil temperature') 
Fitting_type = c(Fitting_type,'Linear fitting')
Data_type = c(Data_type,prefix)

# Quadratic polynomial
mod = lm(Richness~poly(Temp.,2),plot_data) 
AIC = c(AIC,AIC(mod)) %>% round(3)
R2 = c(R2,rsq(mod,adj = T)) %>% round(3)
P_value = c(P_value,lmp(mod)) %>% round(3)
Variable = c(Variable,'Soil temperature')
Fitting_type = c(Fitting_type,'Quadratic polynomial fitting')
Data_type = c(Data_type,prefix)

### pH
# linear
plot_data$Richness = specnumber(t(OTU))
mod = lm(Richness~pH,plot_data) 
AIC = c(AIC,AIC(mod)) %>% round(3)
R2 = c(R2,rsq(mod,adj = T)) %>% round(3)
P_value = c(P_value,lmp(mod)) %>% round(3)
Variable = c(Variable,'pH')
Fitting_type = c(Fitting_type,'Linear fitting')
Data_type = c(Data_type,prefix)

# Quadratic polynomial
mod = lm(Richness~poly(pH,2),plot_data) 
AIC = c(AIC,AIC(mod)) %>% round(3)
R2 = c(R2,rsq(mod,adj = T)) %>% round(3)
P_value = c(P_value,lmp(mod)) %>% round(3)
Variable = c(Variable,'pH')
Fitting_type = c(Fitting_type,'Quadratic polynomial fitting')
Data_type = c(Data_type,prefix)

### Moisture
# linear
plot_data$Richness = specnumber(t(OTU))
mod = lm(Richness~Water_content,plot_data) 
AIC = c(AIC,AIC(mod)) %>% round(3)
R2 = c(R2,rsq(mod,adj = T)) %>% round(3)
P_value = c(P_value,lmp(mod)) %>% round(3)
Variable = c(Variable,'Soil moisture')
Fitting_type = c(Fitting_type,'Linear fitting')
Data_type = c(Data_type,prefix)

# Quadratic polynomial
mod = lm(Richness~poly(Water_content,2),plot_data) 
AIC = c(AIC,AIC(mod)) %>% round(3)
R2 = c(R2,rsq(mod,adj = T)) %>% round(3)
P_value = c(P_value,lmp(mod)) %>% round(3)
Variable = c(Variable,'Soil moisture')
Fitting_type = c(Fitting_type,'Quadratic polynomial fitting')
Data_type = c(Data_type,prefix)

# 99% Identity
load('Project_16S_99')
project = get(paste0('Project_16S_99')) %>% rarefy_even_depth(sample.size=40000)
OTU <- otu_table(project) %>% data.frame() %>% .[,plot_data$Sample_ID]
prefix = '16S_99'

### Temp.
# linear
plot_data$Richness = specnumber(t(OTU))
mod = lm(Richness~Temp.,plot_data) 
AIC = c(AIC,AIC(mod)) %>% round(3)
R2 = c(R2,rsq(mod,adj = T)) %>% round(3)
P_value = c(P_value,lmp(mod)) %>% round(3)
Variable = c(Variable,'Soil temperature')
Fitting_type = c(Fitting_type,'Linear fitting')
Data_type = c(Data_type,prefix)

# Quadratic polynomial
mod = lm(Richness~poly(Temp.,2),plot_data) 
AIC = c(AIC,AIC(mod)) %>% round(3)
R2 = c(R2,rsq(mod,adj = T)) %>% round(3)
P_value = c(P_value,lmp(mod)) %>% round(3)
Variable = c(Variable,'Soil temperature')
Fitting_type = c(Fitting_type,'Quadratic polynomial fitting')
Data_type = c(Data_type,prefix)

### pH
# linear
plot_data$Richness = specnumber(t(OTU))
mod = lm(Richness~pH,plot_data) 
AIC = c(AIC,AIC(mod)) %>% round(3)
R2 = c(R2,rsq(mod,adj = T)) %>% round(3)
P_value = c(P_value,lmp(mod)) %>% round(3)
Variable = c(Variable,'pH')
Fitting_type = c(Fitting_type,'Linear fitting')
Data_type = c(Data_type,prefix)

# Quadratic polynomial
mod = lm(Richness~poly(pH,2),plot_data) 
AIC = c(AIC,AIC(mod)) %>% round(3)
R2 = c(R2,rsq(mod,adj = T)) %>% round(3)
P_value = c(P_value,lmp(mod)) %>% round(3)
Variable = c(Variable,'pH')
Fitting_type = c(Fitting_type,'Quadratic polynomial fitting')
Data_type = c(Data_type,prefix)

### Moisture
# linear
plot_data$Richness = specnumber(t(OTU))
mod = lm(Richness~Water_content,plot_data) 
AIC = c(AIC,AIC(mod)) %>% round(3)
R2 = c(R2,rsq(mod,adj = T)) %>% round(3)
P_value = c(P_value,lmp(mod)) %>% round(3)
Variable = c(Variable,'Soil moisture')
Fitting_type = c(Fitting_type,'Linear fitting')
Data_type = c(Data_type,prefix)

# Quadratic polynomial
mod = lm(Richness~poly(Water_content,2),plot_data) 
AIC = c(AIC,AIC(mod)) %>% round(3)
R2 = c(R2,rsq(mod,adj = T)) %>% round(3)
P_value = c(P_value,lmp(mod)) %>% round(3)
Variable = c(Variable,'Soil moisture')
Fitting_type = c(Fitting_type,'Quadratic polynomial fitting')
Data_type = c(Data_type,prefix)

# 97% Identity
load('project_16S_OTU')
project = get(paste0('project_16S_OTU'))
OTU <- otu_table(project)
plot_data$Soil.Bac.Chao1 = estimateR(t(OTU))['S.chao1',]
plot_data$Richness = specnumber(t(OTU))

### Temp.
# linear
plot_data$Richness = specnumber(t(OTU))
mod = lm(Richness~Temp.,plot_data) 
AIC = c(AIC,AIC(mod)) %>% round(3)
R2 = c(R2,rsq(mod,adj = T)) %>% round(3)
P_value = c(P_value,lmp(mod)) %>% round(3)
Variable = c(Variable,'Soil temperature')
Fitting_type = c(Fitting_type,'Linear fitting')
Data_type = c(Data_type,prefix)

# Quadratic polynomial
mod = lm(Richness~poly(Temp.,2),plot_data) 
AIC = c(AIC,AIC(mod)) %>% round(3)
R2 = c(R2,rsq(mod,adj = T)) %>% round(3)
P_value = c(P_value,lmp(mod)) %>% round(3)
Variable = c(Variable,'Soil temperature')
Fitting_type = c(Fitting_type,'Quadratic polynomial fitting')
Data_type = c(Data_type,prefix)

### pH
# linear
plot_data$Richness = specnumber(t(OTU))
mod = lm(Richness~pH,plot_data) 
AIC = c(AIC,AIC(mod)) %>% round(3)
R2 = c(R2,rsq(mod,adj = T)) %>% round(3)
P_value = c(P_value,lmp(mod)) %>% round(3)
Variable = c(Variable,'pH')
Fitting_type = c(Fitting_type,'Linear fitting')
Data_type = c(Data_type,prefix)

# Quadratic polynomial
mod = lm(Richness~poly(pH,2),plot_data) 
AIC = c(AIC,AIC(mod)) %>% round(3)
R2 = c(R2,rsq(mod,adj = T)) %>% round(3)
P_value = c(P_value,lmp(mod)) %>% round(3)
Variable = c(Variable,'pH')
Fitting_type = c(Fitting_type,'Quadratic polynomial fitting')
Data_type = c(Data_type,prefix)

### Moisture
# linear
plot_data$Richness = specnumber(t(OTU))
mod = lm(Richness~Water_content,plot_data) 
AIC = c(AIC,AIC(mod)) %>% round(3)
R2 = c(R2,rsq(mod,adj = T)) %>% round(3)
P_value = c(P_value,lmp(mod)) %>% round(3)
Variable = c(Variable,'Soil moisture')
Fitting_type = c(Fitting_type,'Linear fitting')
Data_type = c(Data_type,prefix)

# Quadratic polynomial
mod = lm(Richness~poly(Water_content,2),plot_data) 
AIC = c(AIC,AIC(mod)) %>% round(3)
R2 = c(R2,rsq(mod,adj = T)) %>% round(3)
P_value = c(P_value,lmp(mod)) %>% round(3)
Variable = c(Variable,'Soil moisture')
Fitting_type = c(Fitting_type,'Quadratic polynomial fitting')
Data_type = c(Data_type,prefix)

result = data.frame(Data_type,Variable,Fitting_type,R2,P_value,AIC)
write.csv(result,file = paste0(folder,'/Fitting_types_stat.csv'))
