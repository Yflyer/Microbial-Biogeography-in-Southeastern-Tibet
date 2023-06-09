source('0_function.R')
library(lavaan)
library(semTools)
####
#################### save dir #######################
folder = paste0('9_sem')
dir.create(folder)
#################### load data #######################
match_dt <- read.csv('1_Data_Merge/Merge_dt.csv',row.names = 1,stringsAsFactors = T)

#################### bacteria community ###############
st_dt = vroom('7_icamp_ITS_ASV/ITS_ASV.ProcessImportance_EachTurnover.csv')

### stochasticity of each sample
Group = match_dt %>% select(Sample_ID,Plot) 
st1.wg = st_dt %>% 
  left_join(.,Group,by=c('samp1' = 'Sample_ID'))  %>%
  left_join(.,Group,by=c('samp2' = 'Sample_ID'))  %>%
  filter(Plot.x == Plot.y) %>%
  mutate(Sample_ID=samp1)
st2.wg = st_dt %>% 
  left_join(.,Group,by=c('samp1' = 'Sample_ID'))  %>%
  left_join(.,Group,by=c('samp2' = 'Sample_ID'))  %>%
  filter(Plot.x == Plot.y) %>%
  mutate(Sample_ID=samp2)

model_dt = rbind(st1.wg,st2.wg) %>% 
  group_by(Sample_ID) %>% 
  summarise(HeS = mean(HeS),
            HoS = mean(HoS),
            DL = mean(DL),
            HD = mean(HD),
            DR = mean(DR)) %>%
  full_join(.,match_dt) %>%
  select_if(is.numeric) %>% 
  scale() %>%
  as.data.frame()

colnames(model_dt)
############################
# stochasticity and environment change beta effect on Rs patterns
### DL
full_model <- ' 
  Fug_disper ~ pH + SOC + NH4 + NO3 + TN + Temp. + Water_content
  DL ~ Fug_disper + pH + SOC + NH4 + NO3 + TN + Temp. + Water_content
  HD ~ Fug_disper + pH + SOC + NH4 + NO3 + TN + Temp. + Water_content 
  DR ~ Fug_disper + pH + SOC + NH4 + NO3 + TN + Temp. + Water_content 
  Rs ~ pH + SOC + NH4 + NO3 + TN + Temp. + Water_content + Fug_disper + DL + HD + DR
'
# select processes by SEM performance
fit1 <- sem(full_model, data = model_dt)
result <- summary(fit1, standardized = TRUE,fit.measures = TRUE) 
result$PE %>% filter(lhs %in% c('Rs','Fug_disper','DL') & rhs %in% c('DL','Fug_disper'))
fitmeasures(fit1)[c('rmsea','cfi','tli')]

###############################
lm(Rs ~ pH + Fug_disper + SOC + NH4 + NO3 + TN + Temp. + Water_content + DL + HD + DR,model_dt) %>% step(.) %>% summary(.)
lm(DL ~ Fug_disper + pH + SOC + NH4 + NO3 + TN + Temp. + Water_content,model_dt) %>% step(.) %>% summary(.)
lm(Fug_disper ~ pH + SOC + NH4 + NO3 + TN + Temp. + Water_content,model_dt) %>% step(.) %>% summary(.)

# prune paths by stepwise regression
###############################
prune_model <- '
  Fug_disper ~ SOC + Temp.
  DL ~ NO3 + Temp. + Fug_disper
  Rs ~ Fug_disper + TN + Temp. + DL
'
fit1 <- sem(prune_model, data = model_dt)
resid(fit1, "cor")
modificationindices(fit1, minimum.value = 20)
fitmeasures(fit1)[c('rmsea','cfi','tli')]
summary(fit1, standardized = TRUE,fit.measures = TRUE) -> result
write.csv(result$PE,paste0(folder,'/Bac_SEM_pruning.csv'),row.names = F)
write.csv(data.frame(result$FIT),paste0(folder,'/Bac_SEM_pruning_fit.csv'))

