source('0_function.R')
#################### save dir #######################
folder = paste0('S06')
dir.create(folder)

#################### load data ######################
Env = read.csv('S01/MTE_ENV.csv',row.names = 1)
load('project_16S_ASV')
project = get(paste0('project_',prefix))
OTU = otu_table(project) %>% data.frame()
Samp =  group_sum(OTU,sample_data(project)$Plot)
Kgroup = rowSums(Samp>0)

project = prune_taxa(taxa_names(project)[Kgroup>=3], project)
OTU = otu_table(project) %>% data.frame()

All.var = c('Temp.','pH','TN','SOC','Water_content','NH4')

################### loop ###########################
Final_dt = data.frame()
for (j in 1:length(All.var)) {
  r_square=c()
  p_value = c()
  slope = c()
  taxa_ID = c()
  taxaInfo = c()
  for (i in 1:nrow(OTU)) {
    taxa.abund = OTU[i,] %>% as.numeric(.)  #%>% scale(.)
    Fun_value = Env[,All.var[j]] %>% as.numeric(.) #%>% scale(.)
    # match
    Fun_value  = Fun_value[taxa.abund>0]
    taxa.abund = taxa.abund[taxa.abund>0]
    ####
    Cor = cor.test(taxa.abund,Fun_value,method = 'spearman',exact = F)
    
    #r_square=c(r_square,result$adj.r.squared)
    slope = c(slope,Cor$estimate)
    p_value =c(p_value,Cor$p.value)
    taxa_ID = c(taxa_ID,rownames(OTU)[i])
    #taxa_info = c(taxa_info,TAX[rownames(OTU)[i],'Lineage'])
  }
  Result_lm = data.frame(taxa_ID,slope,p_value) # ,taxa_info
  # p adjust
  Result_lm$p_adjust = p.adjust(Result_lm$p_value, method = "BH") 
  Result_lm$K = Kgroup[rownames(OTU)]
  Result_lm[['Variable']]=All.var[j]
  Final_dt=rbind(Result_lm,Final_dt)
  write.csv(Result_lm,file = paste0(folder,'/Cor_test_',All.var[j],'.csv'))
  
}
Final_dt
write.csv(Result_lm,file = paste0(folder,'/Cor_test_all.csv'))

Plot_dt = Final_dt %>%
  drop_na() %>% 
  group_by(K,Variable) %>%
  summarise(Increased = sum(p_value<0.05&slope>0)/n(),
            Decreased = sum(p_value<0.05&slope<0)/n(),
            Notrend = sum(p_value>0.05)/n()) %>%
  gather(key = 'Type',value = 'Proportion',-K,-Variable) %>%
  mutate(K = factor(K,levels=1:12),
         Type = factor(Type,levels=c('Decreased','Notrend','Increased')),
         Variable=factor(Variable,levels=,c('Temp.','pH','Water_content','SOC','TN','NH4')))
write.csv(Plot_dt,file = paste0(folder,'/Cor_test_proportion.csv'))

ggplot(Plot_dt)+geom_col(aes(x=K,y=Proportion,group=Type,fill=Type))+facet_wrap(~Variable,ncol = 3)+
  scale_y_continuous(labels = scales::percent,limits =c(0,1),breaks = seq(0,1,0.25))+
  theme_light()+
  scale_fill_npg()+
  theme(panel.grid = element_blank())+
  force_panelsizes(rows = unit(5, "cm"),cols = unit(3.5, "cm"))
ggsave(filename = paste0(folder,'/Cor_test_proportion.pdf'),dpi=300,width=16,height=14,units='cm')

# K12 temp incr: 26.3 %
# K12 temp decr: 10.9 %