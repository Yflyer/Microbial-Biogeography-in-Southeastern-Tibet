source('0_function.R')
library(ggplot2)
library(ggh4x)
# library(spaa)
#################### save dir #######################
folder = paste0('S05')
dir.create(folder)
########### data ##############
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
# load data
Env = read.csv('S02/MTE_ENV.csv',row.names = 1,stringsAsFactors = T)
Process_bin = read.csv('S01/icamp_16S_ASV/16S_ASV.BinContributeToProcess_EachGroup.csv')
Taxa_Bin = read.csv('S01/icamp_16S_ASV/16S_ASV.Taxon_Bin.csv')

load('project_16S_ASV')
project = get(paste0('project_',prefix))
OTU = otu_table(project) %>% data.frame() %>% .[Taxa_Bin$ID,]
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
# K data
samp =  group_count(OTU,sample_data(project)$Plot)
K = rowSums(samp>0)
Occurrence = rowSums(OTU>0)

# stochastic data
most_common <- function(x) {
  tt <- table(x)
  names(tt[tt == max(tt)])
}

Sto_dt = Taxa_Bin %>% mutate(K, Occurrence) %>%
  group_by(Bin) %>%
  summarise(K = mean(K),
            Occurrence = mean(Occurrence),
            Relative.abundance = sum(TaxonRelativeAbundance),
            Phylum = most_common(phylum),
            Phylum_number = sort(table(phylum), decreasing = TRUE)[1],
            Class = most_common(class),
            Class_number = sort(table(class), decreasing = TRUE)[1],
            Richness = n()) %>%
  mutate(Phylum=gsub('_','',Phylum),
         Class=gsub('_','',Class)) %>%
  filter(Phylum !='')

# Bin process
colnames(Process_bin)
Inner_dt = Process_bin %>%
  filter(Group %in% unique(Env$Elevation)) %>%
  select(-GroupBasedOn,-Method) %>%
  gather(key = Bin,value=ST,-Group,-Process) %>%
  group_by(Bin) %>%
  mutate(Total_ST = sum(ST),
         Bin=gsub('b','B',Bin),
         Type = 'Inner_of_elevation') %>%
  group_by(Bin,Process) %>%
  summarise(ST = sum(ST)/unique(Total_ST))

Between_dt = Process_bin %>%
  filter(!(Group %in% unique(Env$Elevation))) %>%
  select(-GroupBasedOn,-Method) %>%
  gather(key = Bin,value=ST,-Group,-Process) %>%
  group_by(Bin) %>%
  mutate(Total_ST = sum(ST),
         Bin=gsub('b','B',Bin),
         Type = 'Between_elevation') %>%
  group_by(Bin,Process) %>%
  summarise(ST = sum(ST)/unique(Total_ST))

# Merge
plot_dt = Sto_dt %>%
  left_join(Between_dt) %>%
  drop_na()

Index = plot_dt %>% 
  filter(Phylum!='UNKNOWN',Relative.abundance>0) %>%
  group_by(Phylum) %>%
  summarise(Total.abund = sum(Relative.abundance),
            Log.abund = log10(Total.abund))
Order.taxa = Index$Phylum[order(Index$Total.abund,decreasing = T)] %>% .[1:10]
plot_dt$Phylum[!plot_dt$Phylum %in% Order.taxa]='Others'
plot_dt$Phylum = factor(plot_dt$Phylum,levels=c(Order.taxa,'Others'))

### process plot
plot_dt1 = plot_dt %>% filter(Process %in% c('DL'))

stat_dt = plot_dt1 %>% group_by(Process) %>% 
  summarise(Trend = format(lm_par(ST,K)[1], digits = 3),
            Intercept = format(lm_par(ST,K)[3], digits = 3),
            R2 =format(lm_par(ST,K)[4],digits =3,nsmall=3),
            P = labelpstar(lm_par(ST,K)[2])) %>%
  mutate(Text=paste0(Process,' = ',Trend,' * K + ',Intercept,'\nR2 = ',R2,P))
ggplot(plot_dt1)+geom_point(aes(K,ST,color=Phylum))+
  geom_smooth(aes(K,ST),alpha=0.5,method = 'lm',color='black',linewidth=1.5)+
  annotate("text", x = 3, y = 1.1, label = stat_dt$Text,hjust=0,vjust=1)+
  labs(y = "Dispersal limitation(%)")+
  scale_y_continuous(labels = scales::percent,limits = c(0,1.1),breaks = seq(0,1,0.25))+
  theme_bw()+
  guides(fill = 'none')+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        legend.position = 'right',
        legend.key.size = unit(0.5, 'cm'),
        legend.background = element_rect(colour = "black"),
        legend.text = element_text(size=12))+
  force_panelsizes(rows = unit(6, "cm"),cols = unit(6, "cm"))
ggsave(filename = paste0(folder,'/ST_DL.pdf'),dpi=600,width=12,height=12,units='cm')

### process plot
plot_dt1 = plot_dt %>% filter(Process %in% c('HoS'))

stat_dt = plot_dt1 %>% group_by(Process) %>% 
  summarise(Trend = format(lm_par(ST,K)[1], digits = 3),
            Intercept = format(lm_par(ST,K)[3], digits = 3),
            R2 =format(lm_par(ST,K)[4],digits =3,nsmall=3),
            sd =format(lm_par(ST,K)[5],digits =3,nsmall=3),
            P = labelpstar(lm_par(ST,K)[2])) %>%
  mutate(Text=paste0(Process,' = ',Trend,' * K + ',Intercept,'\nR2 = ',R2,P,' ',sd))
ggplot(plot_dt1)+geom_point(aes(K,ST,color=Phylum))+
  geom_smooth(aes(K,ST),alpha=0.5,method = 'lm',color='black',linewidth=1.5)+
  annotate("text", x = 3, y = 1.1, label = stat_dt$Text,hjust=0,vjust=1)+
  labs(y = "Homogenous selection(%)")+
  scale_y_continuous(labels = scales::percent,limits = c(0,1.1),breaks = seq(0,1,0.25))+
  theme_bw()+
  guides(fill = 'none')+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        legend.position = 'right',
        legend.key.size = unit(0.5, 'cm'),
        legend.background = element_rect(colour = "black"),
        legend.text = element_text(size=12))+
  force_panelsizes(rows = unit(6, "cm"),cols = unit(6, "cm"))
ggsave(filename = paste0(folder,'/ST_HoS.pdf'),dpi=600,width=12,height=12,units='cm')

### process plot
plot_dt1 = plot_dt %>% filter(Process %in% c('HeS'))

stat_dt = plot_dt1 %>% group_by(Process) %>% 
  summarise(Trend = format(lm_par(ST,K)[1], digits = 3),
            Intercept = format(lm_par(ST,K)[3], digits = 3),
            R2 =format(lm_par(ST,K)[4],digits =3,nsmall=3),
            P = labelpstar(lm_par(ST,K)[2])) %>%
  mutate(Text=paste0(Process,' = ',Trend,' * K + ',Intercept,'\nR2 = ',R2,P))
ggplot(plot_dt1)+geom_point(aes(K,ST,color=Phylum))+
  geom_smooth(aes(K,ST),alpha=0.5,method = 'lm',color='black',linewidth=1.5)+
  annotate("text", x = 3, y = 0.5, label = stat_dt$Text,hjust=0,vjust=1)+
  labs(y = "Heterogenous selection(%)")+
  scale_y_continuous(labels = scales::percent,limits = c(0,0.5),breaks = seq(0,0.5,0.1))+
  theme_bw()+
  guides(fill = 'none')+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        legend.position = 'right',
        legend.key.size = unit(0.5, 'cm'),
        legend.background = element_rect(colour = "black"),
        legend.text = element_text(size=12))+
  force_panelsizes(rows = unit(6, "cm"),cols = unit(6, "cm"))
ggsave(filename = paste0(folder,'/ST_HeS.pdf'),dpi=600,width=12,height=12,units='cm')

### process plot
plot_dt1 = plot_dt %>% filter(Process %in% c('DR'))

stat_dt = plot_dt1 %>% group_by(Process) %>% 
  summarise(Trend = format(lm_par(ST,K)[1], digits = 3),
            Intercept = format(lm_par(ST,K)[3], digits = 3),
            R2 =format(lm_par(ST,K)[4],digits =3,nsmall=3),
            P = labelpstar(lm_par(ST,K)[2])) %>%
  mutate(Text=paste0(Process,' = ',Trend,' * K + ',Intercept,'\nR2 = ',R2,P))
ggplot(plot_dt1)+geom_point(aes(K,ST,color=Phylum))+
  geom_smooth(aes(K,ST),alpha=0.5,method = 'lm',color='black',linewidth=1.5)+
  annotate("text", x = 3, y = 0.75, label = stat_dt$Text,hjust=0,vjust=1)+
  labs(y = "Heterogenous selection(%)")+
  scale_y_continuous(labels = scales::percent,limits = c(0,0.75),breaks = seq(0,0.75,0.25))+
  theme_bw()+
  guides(fill = 'none')+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        legend.position = 'right',
        legend.key.size = unit(0.5, 'cm'),
        legend.background = element_rect(colour = "black"),
        legend.text = element_text(size=12))+
  force_panelsizes(rows = unit(6, "cm"),cols = unit(6, "cm"))
ggsave(filename = paste0(folder,'/ST_DR.pdf'),dpi=600,width=12,height=12,units='cm')

### process plot
plot_dt1 = plot_dt %>% filter(Process %in% c('HD'))

stat_dt = plot_dt1 %>% group_by(Process) %>% 
  summarise(Trend = format(lm_par(ST,K)[1], digits = 3),
            Intercept = format(lm_par(ST,K)[3], digits = 3),
            R2 =format(lm_par(ST,K)[4],digits =3,nsmall=3),
            P = labelpstar(lm_par(ST,K)[2])) %>%
  mutate(Text=paste0(Process,' = ',Trend,' * K + ',Intercept,'\nR2 = ',R2,P))
ggplot(plot_dt1)+geom_point(aes(K,ST,color=Phylum))+
  geom_smooth(aes(K,ST),alpha=0.5,method = 'lm',color='black',linewidth=1.5)+
  annotate("text", x = 3, y = 0.2, label = stat_dt$Text,hjust=0,vjust=1)+
  labs(y = "Heterogenous selection(%)")+
  scale_y_continuous(labels = scales::percent,limits = c(0,0.2),breaks = seq(0,0.2,0.05))+
  theme_bw()+
  guides(fill = 'none')+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        legend.position = 'right',
        legend.key.size = unit(0.5, 'cm'),
        legend.background = element_rect(colour = "black"),
        legend.text = element_text(size=12))+
  force_panelsizes(rows = unit(6, "cm"),cols = unit(6, "cm"))
ggsave(filename = paste0(folder,'/ST_HD.pdf'),dpi=600,width=12,height=12,units='cm')

  