source('0_function.R')
#################### save dir #######################
folder = paste0('1_Data_Merge')
dir.create(folder)

load('project_16S_ASV')
project = get(paste0('project_',prefix)) 

OTU <- otu_table(project)
TAX <- tax_table(project)
ENV <- sample_data(project)

# Bacterial alpha
Merge_dt = data.frame(ENV)
Merge_dt$Sample_ID = rownames(ENV)
Merge_dt$Bac_shannon=ENV$Shannon
Merge_dt$Bac_richness=ENV$Richness
Merge_dt$Bac_PD=ENV$PD

# Bacterial beta
Elevation = Merge_dt$Elevation %>% as.factor()
OTU.dist = t(OTU) %>% vegdist()
Merge_dt$Bac_disper = betadisper(OTU.dist,Elevation) %>% .$distance

load('project_ITS_ASV')
project = get(paste0('project_',prefix)) 

OTU <- otu_table(project)
ENV <- sample_data(project)

# Fungal alpha
Merge_dt$Fug_shannon=ENV$Shannon
Merge_dt$Fug_richness=ENV$Richness
Merge_dt$Fug_PD=ENV$PD

# Fungal beta
OTU.dist = t(OTU) %>% vegdist()
Merge_dt$Fug_disper = betadisper(OTU.dist,Elevation) %>% .$distance

# rename table
colnames(Merge_dt)
Merge_dt = Merge_dt %>% select(-Evenness,-Richness,-Invsimpson,-Simpson,-PD,-Shannon)

write.csv(Merge_dt,file = paste0(folder,'/Merge_dt.csv'))
