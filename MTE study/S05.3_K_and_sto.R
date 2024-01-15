source('0_function.R')
library(ggplot2)
library(ggh4x)
# library(spaa)
#################### save dir #######################
folder = paste0('S05')
dir.create(folder)
########### data ##############
# load data
Env = read.csv('S01/MTE_ENV.csv',row.names = 1,stringsAsFactors = T)
Process_bin = read.csv('S01/icamp_16S_ASV/16S_ASV.BinContributeToProcess_EachGroup.csv')
Taxa_Bin = read.csv('S01/icamp_16S_ASV/16S_ASV.Taxon_Bin.csv')

load('project_16S_ASV')
project = get(paste0('project_',prefix))
OTU = otu_table(project) %>% data.frame() %>% .[Taxa_Bin$ID,]

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
  left_join(Inner_dt) %>%
  drop_na()

Index = plot_dt %>% 
  filter(Phylum!='UNKNOWN',Relative.abundance>0) %>%
  group_by(Phylum) %>%
  summarise(Total.abund = sum(Relative.abundance),
            Log.abund = log10(Total.abund))
Order.taxa = Index$Phylum[order(Index$Total.abund,decreasing = T)] %>% .[1:10]

summary(plot_dt$K)
length(plot_dt$K[plot_dt$K<6.576])

plot_dt$Phylum[!plot_dt$Phylum %in% Order.taxa]='Others'
plot_dt$Phylum = factor(plot_dt$Phylum,levels=c(Order.taxa,'Others'))
plot_dt$K_type = 'K>7'
plot_dt$K_type[plot_dt$K<=7]  = 'K<=7'
plot_dt$K_type = factor(plot_dt$K_type,c('K<=7','K>7'))
plot_dt$Process = factor(plot_dt$Process,levels = c('DL','DR','HD','HeS','HoS'))

ggplot(plot_dt)+
  geom_boxplot(aes(x=K_type,y=ST,group=K_type,fill=K_type))+
  #scale_y_continuous(labels = scales::percent,limits = c(0,1.1),breaks = seq(0,1,0.5))+
  facet_wrap(~Process,ncol = 1,scales = 'free')+
  theme_bw()+
  guides(color = 'none')+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        legend.position = 'right',
        legend.key.size = unit(0.5, 'cm'),
        legend.background = element_rect(colour = "black"),
        legend.text = element_text(size=12))+
  force_panelsizes(rows = unit(3, "cm"),cols = unit(4, "cm"))
ggsave(filename = paste0(folder,'/Processs_K_type_difference.pdf'),dpi=600,width=12,height=18,units='cm')

plot_dt$Process_type = 'Stochastic'
plot_dt$Process_type[plot_dt$Process%in%c('HoS','HeS')]  = 'Deterministic'
plot_dt2 = plot_dt
plot_dt2$K_type='The whole community'
total_dt =rbind(plot_dt,plot_dt2) %>% group_by(K_type,Process_type,Bin) %>%
  summarise(ST = sum(ST))
total_dt$K_type = factor(total_dt$K_type,c('K<=7','K>7','The whole community'))

ggplot(total_dt)+
  geom_boxplot(aes(x=Process_type,y=ST,fill=K_type))+
  #scale_y_continuous(labels = scales::percent,limits = c(0,1.1),breaks = seq(0,1,0.5))+
  theme_bw()+
  guides(color = 'none')+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        legend.position = 'right',
        legend.key.size = unit(0.5, 'cm'),
        legend.background = element_rect(colour = "black"),
        legend.text = element_text(size=12))+
  facet_wrap(~Process_type,ncol = 1,scales = 'free')+
  force_panelsizes(rows = unit(7, "cm"),cols = unit(4, "cm"))
ggsave(filename = paste0(folder,'/Overall_Processs_K_type_difference.pdf'),dpi=600,width=12,height=12,units='cm')
