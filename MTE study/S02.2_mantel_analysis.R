source('0_function.R')
#################### save dir #######################
folder = paste0('s02')
dir.create(folder)
###
Dist.list = readRDS(paste0(folder,'/Dist.list'))
##########################################
### stream habitat
plot_data = read.csv('S01/MTE_ENV.csv',row.names = 1)
colnames(plot_data)
### create var
All.var = c('Temp.','pH','TN','SOC','Water_content','NH4')
Result_dt = data.frame()
### bac
for (n in 1:length(Dist.list)) {
  mantel_r = p_value = c()
  for (i in 1:length(All.var)) {
    Env.ctr = plot_data[,All.var[(-i)] ]  %>% scale() %>% dist()
    Env.var = plot_data[,All.var[(i)] ] %>% scale() %>% dist()
    Result = mantel.partial(Dist.list[[n]],Env.var,Env.ctr,method = "spearman")
    mantel_r =c(mantel_r,Result$statistic)
    p_value =c(p_value,Result$signif)
  }
  sub_result = data.frame(Factors=All.var,mantel_r,p_value)
  sub_result$Habitat = str_split(names(Dist.list)[n],'[.]')[[1]][1]
  sub_result$Taxa_group = str_split(names(Dist.list)[n],'[.]')[[1]][2]
  sub_result$Beta = str_split(names(Dist.list)[n],'[.]')[[1]][3]
  Result_dt = rbind(Result_dt,sub_result)
}

write.csv(Result_dt,file = paste0(folder,'/Mantel_result_soil.csv'))
##########################################

