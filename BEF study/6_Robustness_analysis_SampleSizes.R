source('0_function.R')

#################### save dir #######################
folder = paste0('6_Robustness_analysis_SampleSizes/')
dir.create(folder)
############### Env data pre
ENV <- read.csv('1_Data_Merge/Merge_dt.csv',row.names = 1,stringsAsFactors = T)
Elevation = ENV$Elevation %>% as.factor()

########################
load('project_16S_ASV')
project=get(paste0('project_',prefix))
OTU <- otu_table(project)
TAX <- tax_table(project) %>% data.frame()

OTU.dist = t(OTU) %>% vegdist(.)
ENV$Disper = betadisper(OTU.dist,Elevation) %>% .$distance
#########################
Final_result=data.frame()
test_size=3:9

######################## Parallel back-end
library(doParallel) 
cl = makeCluster(detectCores() ) 
## register the cluster 
registerDoParallel(cl) 
for (i in test_size) {
  Boots_pool = combn(1:11,i) 
  ## do parallel computaitons with foreach 
  P_value = foreach(n=1:ncol(Boots_pool),.combine = c,.packages = c("dplyr")) %dopar% { 
    Boots_seq = Boots_pool[,n]
    Dt = ENV %>% group_by(Plot) %>% mutate(Sample_rank=c(1:11)) %>% filter(Sample_rank %in% Boots_seq)
    #### raw effect
    result = lm(Rs~Disper,data=Dt) %>% summary() 
    result$coefficients %>% .['Disper','Pr(>|t|)']
  } 
  Result=data.frame(Sample_size=i,Boots_times=ncol(Boots_pool),P_value)
  Final_result=rbind(Final_result,Result)
}

Final_result$Sign = 'Unsign.'
Final_result$Sign[Final_result$P_value<0.05] = 'Sign.'

N = annotate("text", y=4.5,x=as.factor(3:9),label= paste('N = ',3:9))
ggplot(Final_result)+geom_jitter(aes(x=as.factor(Sample_size),y=-log10(P_value),color=Sign),size=1)+N
ggsave(filename = paste0(folder,'/Bac_boots.pdf'),dpi=900,width=14,height=12,units='cm')
write.csv(Final_result,paste0(folder,'/Bac_boots.csv'))

########################
load('project_ITS_ASV')
project=get(paste0('project_',prefix))
OTU <- otu_table(project)
TAX <- tax_table(project) %>% data.frame()
ENV <- sample_data(project) %>% data.frame()
OTU.dist = t(OTU) %>% vegdist(.)
Elevation=ENV$Elevation
ENV$Disper = betadisper(OTU.dist,Elevation) %>% .$distance
#########################
Final_result=data.frame()
test_size=3:9
for (i in test_size) {
  Boots_pool = combn(1:11,i) 
  ## do parallel computaitons with foreach 
  P_value = foreach(n=1:ncol(Boots_pool),.combine = c,.packages = c("dplyr")) %dopar% { 
    Boots_seq = Boots_pool[,n]
    Dt = ENV %>% group_by(Plot) %>% mutate(Sample_rank=c(1:11)) %>% filter(Sample_rank %in% Boots_seq)
    #### raw effect
    result = lm(Rs~Disper,data=Dt) %>% summary() 
    result$coefficients %>% .['Disper','Pr(>|t|)']
  } 
  Result=data.frame(Sample_size=i,Boots_times=ncol(Boots_pool),P_value)
  Final_result=rbind(Final_result,Result)
}
Final_result$Sign = 'Unsign.'
Final_result$Sign[Final_result$P_value<0.05] = 'Sign.'

N = annotate("text", y=4.5,x=as.factor(3:9),label= paste('N = ',3:9))
ggplot(Final_result)+geom_jitter(aes(x=as.factor(Sample_size),y=-log10(P_value),color=Sign),size=1)+N
ggsave(filename = paste0(folder,'/Fug_boots.pdf'),dpi=900,width=14,height=12,units='cm')
write.csv(Final_result,paste0(folder,'/Fug_boots.csv'))

stopCluster(cl)
