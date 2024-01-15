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
st_dt = vroom('7_icamp_16S_ASV/16S_ASV.ProcessImportance_EachTurnover.csv')

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
  Bac_disper ~ pH + SOC + NH4 + NO3 + TN + Temp. + Water_content
  DL ~ Bac_disper + pH + SOC + NH4 + NO3 + TN + Temp. + Water_content
  HD ~ Bac_disper + pH + SOC + NH4 + NO3 + TN + Temp. + Water_content 
  DR ~ Bac_disper + pH + SOC + NH4 + NO3 + TN + Temp. + Water_content 
  Rs ~ pH + SOC + NH4 + NO3 + TN + Temp. + Water_content  + Bac_disper + DL + HD + DR
'

# select processes by SEM performance
fit1 <- sem(full_model, data = model_dt)
result <- summary(fit1, standardized = TRUE,fit.measures = TRUE) 
result$PE %>% filter(lhs %in% c('Rs','Bac_disper','DL') & rhs %in% c('DL','Bac_disper'))


###############################
lm(Rs ~ pH + SOC + NH4 + NO3 + TN + Temp. + Water_content  + Bac_disper + DL + HD + DR,model_dt) %>% step(.) %>% summary(.)
lm(DL ~ Bac_disper + pH + SOC + NH4 + NO3 + TN + Temp. + Water_content,model_dt) %>% step(.) %>% summary(.)
lm(HD ~ Bac_disper + pH + SOC + NH4 + NO3 + TN + Temp. + Water_content,model_dt) %>% step(.) %>% summary(.)
lm(Bac_disper ~ pH + SOC + NH4 + NO3 + TN + Temp. + Water_content,model_dt) %>% step(.) %>% summary(.)

# prune paths by stepwise regression
###############################
prune_model <- '
  Bac_disper ~ NO3 + Temp. + Water_content
  DL ~ Bac_disper + NH4 + Temp.
  HD ~ Bac_disper + NO3 + Temp.
  Rs ~ DL + HD + TN + Temp. + Water_content
'
fit1 <- sem(prune_model, data = model_dt)
resid(fit1, "cor")
modificationindices(fit1, minimum.value = 20)
fitmeasures(fit1)[c('rmsea','cfi','tli')]
summary(fit1, standardized = TRUE,fit.measures = TRUE) -> result
write.csv(result$PE,paste0(folder,'/Bac_disper_SEM_pruning.csv'),row.names = F)
write.csv(data.frame(result$fit),paste0(folder,'/Bac_disper_SEM_pruning_fit.csv'))

