source('0_function.R')
library(lme4)
library(rsq)
library(car) # mixed model anova
library(ggplot2)
library(ggh4x)
library(viridis)
#################### save dir #######################
folder = paste0('S03')
dir.create(folder)
##################################
Env = read.csv('S01/MTE_ENV.csv',row.names = 1)

### Load phylo signal data
load('project_16S_ASV')
dt = read.csv('S02/Phylo_Kgroup_LocalScale_16S_ASV_layers.csv',row.names = 1,stringsAsFactors = T)
dt = left_join(dt,Env) 

# initialize stat data
dt$Elevation = as.factor(dt$Elevation)
dt$K = as.factor(dt$K)
stat_data = dt %>% group_by(K) %>% mutate_if(is.numeric, scale)

### Richness driver
Model_summary = data.frame()

k_list=unique(stat_data$K)
for (k in k_list) {
  sub_data = filter(stat_data,K==k)
  
  formula.control = 'Richness~HeE+Temp.+pH+SOC+TN+CN_ratio'
  mlm <- lm(formula.control, data = sub_data)
  slm = step(mlm,direction = "backward")
  #vif(slm)
  sum.slm = summary(slm)
  rsq.result = rsq.partial(slm,adj = T)
  
  sub_result = data.frame(Variable = rsq.result$variable, Partial.R2 = rsq.result$partial.rsq, P = sum.slm$coefficients[-1,'Pr(>|t|)'],K=k)
  Model_summary = rbind(Model_summary,sub_result)
  
}
Model_summary$Star = sapply(Model_summary$P,labelpstar)
Model_summary$Sign = 'Significant'
Model_summary$Sign[Model_summary$P>=0.05] = 'Insignificant'
Model_summary$K = factor(Model_summary$K,levels = c(1:12))
Model_summary$Variable = factor(Model_summary$Variable,levels = c('HeE','pH','CN_ratio','SOC','TN','Temp.'))
ggplot(Model_summary)+geom_col(aes(K,Partial.R2,group=Variable,fill=Variable))+
  scale_fill_manual(values = c(hcl.colors(7, "RdYlGn"),'#FFFFFF'))+
  theme_light()+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14))+
  force_panelsizes(rows = unit(6, "cm"),cols = unit(6, "cm"))
ggsave(filename = paste0(folder,'/K_richness_stepwise_R2.pdf'),dpi=600,width=14,height=14,units='cm')

# add mpd and mntd
Model_summary = data.frame()
k_list=unique(stat_data$K)
for (k in k_list) {
  sub_data = filter(stat_data,K==k)
  
  formula.control = 'mpd~HeE+Temp.+pH+SOC+TN+CN_ratio'
  mlm <- lm(formula.control, data = sub_data)
  slm = step(mlm,direction = "backward")
  #vif(slm)
  sum.slm = summary(slm)
  rsq.result = rsq.partial(slm,adj = T)
  
  sub_result = data.frame(Variable = rsq.result$variable, Partial.R2 = rsq.result$partial.rsq, P = sum.slm$coefficients[-1,'Pr(>|t|)'],K=k)
  Model_summary = rbind(Model_summary,sub_result)
  
}
Model_summary$Star = sapply(Model_summary$P,labelpstar)
Model_summary$Sign = 'Significant'
Model_summary$Sign[Model_summary$P>=0.05] = 'Insignificant'
Model_summary$K = factor(Model_summary$K,levels = c(1:12))
Model_summary$Variable = factor(Model_summary$Variable,levels = c('HeE','pH','CN_ratio','SOC','TN','Temp.','mpd','mntd'))
ggplot(Model_summary)+geom_col(aes(K,Partial.R2,group=Variable,fill=Variable))+
  scale_fill_manual(values = c(hcl.colors(7, "RdYlGn"),'grey20','grey40'))+
  theme_light()+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14))+
  force_panelsizes(rows = unit(6, "cm"),cols = unit(6, "cm"))
ggsave(filename = paste0(folder,'/K_mpd_stepwise_R2.pdf'),dpi=600,width=14,height=14,units='cm')

Model_summary = data.frame()
k_list=unique(stat_data$K)
for (k in k_list) {
  sub_data = filter(stat_data,K==k)
  
  formula.control = 'mntd~HeE+Temp.+pH+SOC+TN+CN_ratio'
  mlm <- lm(formula.control, data = sub_data)
  slm = step(mlm,direction = "backward")
  #vif(slm)
  sum.slm = summary(slm)
  rsq.result = rsq.partial(slm,adj = T)
  
  sub_result = data.frame(Variable = rsq.result$variable, Partial.R2 = rsq.result$partial.rsq, P = sum.slm$coefficients[-1,'Pr(>|t|)'],K=k)
  Model_summary = rbind(Model_summary,sub_result)
  
}
Model_summary$Star = sapply(Model_summary$P,labelpstar)
Model_summary$Sign = 'Significant'
Model_summary$Sign[Model_summary$P>=0.05] = 'Insignificant'
Model_summary$K = factor(Model_summary$K,levels = c(1:12))
Model_summary$Variable = factor(Model_summary$Variable,levels = c('HeE','pH','CN_ratio','SOC','TN','Temp.','mpd','mntd'))
ggplot(Model_summary)+geom_col(aes(K,Partial.R2,group=Variable,fill=Variable))+
  scale_fill_manual(values = c(hcl.colors(7, "RdYlGn"),'grey20','grey40'))+
  theme_light()+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14))+
  force_panelsizes(rows = unit(6, "cm"),cols = unit(6, "cm"))
ggsave(filename = paste0(folder,'/K_mntd_stepwise_R2.pdf'),dpi=600,width=14,height=14,units='cm')
