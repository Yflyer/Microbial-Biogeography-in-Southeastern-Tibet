source('0_function.R')
library(ggh4x)
library(rsq)
#################### save dir #######################
folder = paste0('S01')
dir.create(folder)

#################### Sample size effect
lmp = function (modelobject) {
  #if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

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

bootstrap_test = function(i,plot_data,project,Sample_size){
  require(phyloseq)
  require(dplyr)
  require(rsq)
  require(vegan)
  sub_data = plot_data %>% group_by(Elevation) %>% sample_n(Sample_size)
  sub_project =  prune_samples(plot_data$Sample_ID %in% sub_data$Sample_ID,project)
  OTU <- otu_table(sub_project)
  
  # Sort
  sub_data = sub_data %>% arrange(factor(Sample_ID, levels = colnames(OTU)))
  
  sub_data$Soil.Bac.Chao1 = estimateR(t(OTU))['S.chao1',]
  sub_data$Soil.Bac.Richness = specnumber(t(OTU))
  # Ea
  Fit_par = lm_par(log(sub_data$Soil.Bac.Chao1),sub_data$Inv.kT)
  Ea = Fit_par[1] %>% round(.,3)
  # Temp
  Fit_par = lm_par(sub_data$Soil.Bac.Richness,sub_data$Temp.)
  Eestimate.temp = Fit_par[1] %>% round(.,3)
  Eestimate.sd.temp=Fit_par[5] %>% round(.,3)
  P.temp = Fit_par[2] #%>% round(.,3)
  R2.temp = Fit_par[4] %>% round(.,3)
  # Ele
  Fit_par = lm_par(sub_data$Soil.Bac.Richness,as.numeric(sub_data$Elevation))
  Eestimate.ele = Fit_par[1] %>% round(.,3)
  Eestimate.sd.ele=Fit_par[5] %>% round(.,3)
  P.ele = Fit_par[2] #%>% round(.,3)
  R2.ele = Fit_par[4] %>% round(.,3)
  # # multi
  # formula.control = paste0('Soil.Bac.Richness~HeE+Temp.+Water_content+pH+SOC+TN+CN_ratio+NH4') %>% as.formula()
  # mlm <- lm(formula.control, data = sub_data)
  # slm = step(mlm,direction = "backward")
  # sum.slm = summary(slm)
  # rsq.result = rsq.partial(slm,adj = T)
  # P.R2.prop = (rsq.result$partial.rsq[2] / sum.slm$adj.r.squared) %>% round(.,3)
  
  Log = data.frame(Sample_size = Sample_size,
                   Bootstrap_num = i, 
                   #Replicate_log = n,
                   Sample_list = paste(sub_data$Sample_ID,collapse = ' '),
                   Eestimate.ele,Eestimate.sd.ele,P.ele,R2.ele,
                   Eestimate.temp,Eestimate.sd.temp,P.temp,R2.temp
                   # P.R2.prop,Ea
                   )
  Log
}
####################
load('project_16S_ASV')
project = get(paste0('project_16S_ASV'))
plot_data = read.csv('S01/MTE_ENV.csv',row.names = 1) #%>% filter(Elevation!=2013)

Elevation=plot_data$Elevation
plot_data$Elevation = as.factor(plot_data$Elevation)

dt = data.frame()
size_list = c(3:10)
Final_dt = data.frame()

# parallel test
library(snowfall)
# library(data.table)
sfInit(parallel=T,cpus = 12)
sfExport('bootstrap_test','lm_par','lmp','plot_data','project')

for (Sample_size in size_list) {
  Log = data.frame()
  random_list = list()
  sfExport('Sample_size')
  Result = sfLapply(1:100,function(i)bootstrap_test(i,plot_data,project,Sample_size))
  Result = Reduce(rbind,Result) %>% data.frame()
  Final_dt = rbind(Final_dt,Result)
}

# library(data.table)
# Final_dt = rbindlist(Result)
# 
#colnames(Final_dt)= c('Sample_size' ,'Bootstrap_num', 'Sample_list', 'Eestimate.ele','Eestimate.sd.ele','P.ele','R2.ele', 'Eestimate.temp','Eestimate.sd.temp','P.temp','R2.temp')
Final_dt = Final_dt %>% distinct(Sample_list, .keep_all = T)
Final_dt$Sign = 'Unsign.'
Final_dt$Sign[Final_dt$P.temp<0.05] = 'Sign.'
Final_dt$P.temp = unlist(Final_dt$P.temp)
N = annotate("text", y=15,x=as.factor(3:10),label= paste('N = ',3:10))
ggplot(Final_dt)+geom_jitter(aes(x=as.factor(Sample_size),y=-log10(P.temp),color=Sign),size=1)+N
ggsave(filename = paste0(folder,'/MTE_bootstrap_examinations.pdf'),dpi=900,width=14,height=12,units='cm')

write.csv(Final_dt,file = paste0(folder,'/MTE_bootstrap_examinations.csv'))

