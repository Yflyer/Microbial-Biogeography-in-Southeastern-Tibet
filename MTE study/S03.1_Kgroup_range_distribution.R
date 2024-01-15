source('0_function.R')
#################### save dir #######################
folder = paste0('S03')
dir.create(folder)
########### input ##############
# set input
Env = read.csv('S01/MTE_ENV.csv',row.names = 1)
load('project_16S_ASV')
project = get(paste0('project_',prefix))
OTU = otu_table(project) %>% data.frame()
Samp =  group_sum(OTU,sample_data(project)$Plot)
Kgroup = rowSums(Samp>0)

library(snowfall)
sfInit(cpus = 8,parallel = 8)
sfExport("Samp","OTU","Kgroup","Env")
sfLibrary(dplyr)

get_rangedt = function(i,OTU,Samp,Kgroup,Env){
  Site_list = colnames(OTU)[OTU[i,]>0]
  Max.abund.site = which.max(OTU[i,]) %>% names()
  dt = Env %>% filter(Sample_ID %in% Site_list) %>%
    summarise(OTU_ID = rownames(OTU)[i],
              K = Kgroup[i],
              Max.abund = max(OTU[i,]),
              Max.abund.ele = max(Samp[i,]),
              Mean.abund.ele = mean(Samp[i,][Samp[i,]>0]),
              Mid.ele = (max(Elevation)+min(Elevation))/2,
              Max.ele = max(Elevation),
              Min.ele = min(Elevation),
              Range.ele = max(Elevation)-min(Elevation),
              Abund.ele = Elevation[Sample_ID==Max.abund.site],
              Mid.temp = (max(Temp.)+min(Temp.))/2,
              Max.temp = max(Temp.),
              Min.temp = min(Temp.),
              Range.temp = max(Temp.)-min(Temp.),
              Abund.temp = Temp.[Sample_ID==Max.abund.site]
              )
  return(dt)
}

Krange_dt = sfLapply(1:nrow(OTU),get_rangedt,OTU,Samp,Kgroup,Env)
sfStop()

write.csv(Krange_dt,file = paste0(folder,'/Krange_distribution.csv'))

# make plot of temprange distribution
library(data.table)
plot_dt = rbindlist(Krange_dt) %>% 
  filter(K>=3) %>%
  group_by(K) %>%
  mutate(Richness_K = n()) %>%
  group_by(K,Abund.ele) %>%
  summarise(Percentage = n()/unique(Richness_K),
            Mid.temp = mean(Mid.temp),
            Min.temp = mean(Min.temp),
            Max.temp = mean(Max.temp),
            Range.temp = mean(Range.temp),
            Mean.abund.ele=mean(Mean.abund.ele))

plot_dt$K = as.factor(plot_dt$K)
ggplot(plot_dt, aes(x = Abund.ele, Mid.temp)) +
  #geom_errorbar(aes(ymin = Min.temp, ymax = Max.temp,color=K), width = 0.1)+
  geom_point(aes(fill=K,size = Percentage),shape=21,alpha=0.8) +
  #geom_smooth(method = "lm", se = FALSE)+
  theme_bw()
ggsave(filename = paste0(folder,'/Krange_distribution.pdf'),dpi=600,width=18,height=14,units='cm')


ggplot(plot_dt, aes(x = Abund.ele, Range.temp)) +
  #geom_errorbar(aes(ymin = Min.temp, ymax = Max.temp,color=K), width = 0.1)+
  geom_point(aes(fill=K,size = Percentage),shape=21,alpha=0.8) +
  #geom_smooth(method = "lm", se = FALSE)+
  theme_bw()

# pie plot of K group percentage
plot_dt = rbindlist(Krange_dt) %>% 
  group_by(K) %>%
  summarise(Percentage = n()/nrow(OTU)) %>%
  filter(K>=3)
plot_dt = rbind(plot_dt,c(NA,1-sum(plot_dt$Percentage)))
plot_dt$K = as.factor(plot_dt$K)

ggplot(plot_dt, aes(x = "", y = Percentage, fill = K)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_void()+
  geom_text(aes(label = scales::percent(Percentage)), position = position_stack(vjust = 0.5))
ggsave(filename = paste0(folder,'/Kgroup_percentage.pdf'),dpi=600,width=18,height=14,units='cm')

