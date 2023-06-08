# AO-plot
library(ggplot2)
library(ggh4x)
source('0_function.R')
#################### save dir #######################
folder = paste0('8_A&O_analysis')
dir.create(folder)
############### Env data pre
ENV <- read.csv('1_Data_Merge/Merge_dt.csv',row.names = 1,stringsAsFactors = T)
Elevation = ENV$Elevation %>% as.factor()

################### Bac ############################
### load project  
load('project_16S_ASV')
project = get(paste0('project_',prefix)) 

OTU <- otu_table(project)
#########################################################################
####count occupancy frequency on sample and plot levels
count_table=group_count(table = OTU, group = Elevation) %>% as.data.frame(.)
Occurrence = rowSums(count_table)
Abundance=rowSums(OTU)/rowSums(OTU>0)
Species_occupancy=factor(rowSums(count_table>0))

# get plot occurrence according to specific frequency condition
ab_oc_data=data.frame(Occurrence,Abundance,Species_occupancy,Prefix='Bacterial community')

################### Fug ############################
### load project  
load('project_ITS_ASV')
project = get(paste0('project_',prefix)) 
OTU <- otu_table(project)
#########################################################################
####count occupancy frequency on sample and plot levels
count_table=group_count(table = OTU, group = Elevation) %>% as.data.frame(.)
Occurrence = rowSums(count_table)
Abundance=rowSums(OTU)/rowSums(OTU>0)
Species_occupancy=factor(rowSums(count_table>0))

# get plot occurrence according to specific frequency condition
ab_oc_data=rbind(ab_oc_data,data.frame(Occurrence,Abundance,Species_occupancy,Prefix='Fungal community'))
# data transformation 
# optional if you used reads data but not percentage)
ab_oc_data$Log.abundance = log1p(ab_oc_data$Abundance) 
ab_oc_data$Taxon_ID = factor(rownames(ab_oc_data))
ab_oc_data$Species_occupancy = factor(ab_oc_data$Species_occupancy,levels = 1:12)

################################################################
sum = ab_oc_data %>% 
  group_by(Species_occupancy,Prefix) %>%
  summarise(Avg.abundance = mean(Abundance),
            sd.abundance = sd(Abundance),
            Species.number = length(Species_occupancy))

#########################################################################
PTs <- theme(legend.title = element_blank(),
             legend.position="bottom",
             axis.text = element_text(size=12),
             panel.background=element_blank(),
             panel.border = element_rect(colour = "grey", fill=NA, linewidth=1))
##############################
sfm = scale_fill_manual(values = colorRampPalette(colors = c('#11659a','#d0dfe6','#d8e3e7','#f1c4cd','#ec7696'))(12))
scm = scale_color_manual(values = colorRampPalette(colors = c('#11659a','#d0dfe6','#d8e3e7','#f1c4cd','#ec7696'))(12))
ssm = scale_shape_manual(values = c(16,21))

plot_data = ab_oc_data %>% 
  group_by(Prefix,Species_occupancy) %>%
  filter(!is.outlier(Log.abundance, thres = 2)) # this step only improves the visualization. Do not use for statistic

# AO pattern
ggplot(plot_data)+
  geom_boxplot(aes(Species_occupancy,Log.abundance,fill=Species_occupancy),outlier.shape = NA)+
  sfm+
  labs(y = "Log(Abundance)", x = "Occupancy")+
  #scale_y_continuous(limits = c(0,5),breaks = seq(0,5,2))+
  facet_grid(Prefix~.,scales = 'free_y') +
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.spacing = unit(0.2, "lines"),
        legend.key.size = unit(0.1, "cm"),
        legend.key.width = unit(0.5, "cm"))+
  force_panelsizes(rows = unit(2.5, "cm"),cols = unit(6.5, "cm"))
ggsave(filename = paste0(folder,'/AO_abundance.pdf'),dpi=300,width=14,height=8,units='cm')
