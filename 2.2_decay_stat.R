source('0_function.R')
library(rsq)
########################
folder = paste0('2_decay')
dir.create(folder)
########################
decay_data = read.csv('2_decay/decay_data.csv')

########################
### environmental factors
colnames(decay_data)
Spatialfactor = c('Temp.','pH')
Nutrient = c('Water_content','SOC','TN','NO3','NH4')
Var=c(Spatialfactor,Nutrient)

######################## regional scale
Scale = 'Regional'

### decay estimation
dt = decay_data %>% select_if(is.numeric) %>% scale() %>% as.data.frame()

################# Bac
# distance only
decay <- lm(Bac_Bray~Distance, data=dt)
summary(decay) 

# full model
Formula =paste(c('Bac_Bray~Distance',Var),collapse='+')
decay <- lm(Formula, data=dt)
#
anova(decay)
partial.r2 = rsq.partial(decay,adj = T)
Model_summary = summary(decay) 
result = Model_summary$coefficients %>% .[-1,] %>% as.data.frame()
Cof_table = data.frame(Variable = rownames(result),
                       Coefficient = round(result$Estimate,3), 
                       Sign = sapply(result$`Pr(>|t|)`,labelpstar) ,
                       partial.r2 = round(partial.r2$partial.rsq,3))
write.csv(Cof_table,file = paste0(folder,'/Cof_decay_bac.csv'))

############## fug
# distance only
decay <- lm(Bac_Bray~Distance, data=dt)
summary(decay) 

# full_model
Formula =paste(c('Fug_Bray~Distance',Var),collapse='+')
decay <- lm(Formula, data=dt)
#
anova(decay)
partial.r2 = rsq.partial(decay,adj = F)
Model_summary = summary(decay) 
result = Model_summary$coefficients %>% .[-1,] %>% as.data.frame()
Cof_table = data.frame(Variable = rownames(result),
                       Coefficient = round(result$Estimate,3), 
                       Sign = sapply(result$`Pr(>|t|)`,labelpstar) ,
                       partial.r2 = round(partial.r2$partial.rsq,3))
write.csv(Cof_table,file = paste0(folder,'/Cof_decay_fug.csv'))
