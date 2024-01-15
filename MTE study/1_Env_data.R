source('0_function.R')
#################### save dir #######################
folder = paste0('1_Env_and_Comm')
dir.create(folder)

load('project_16S_ASV')
project = get(paste0('project_',prefix)) 
project = subset_samples(project, Elevation !=2013)

OTU <- otu_table(project)
TAX <- tax_table(project)
ENV <- sample_data(project)

ENV_dt = data.frame(ENV)
ENV_dt$label = rownames(ENV)
ENV_dt$Bac_shannon=ENV$Shannon
ENV_dt$Bac_invsimpson=ENV$Invsimpson
ENV_dt$Bac_richness=ENV$Richness
ENV_dt$Bac_PD=ENV$PD

load('project_ITS_ASV')
project = get(paste0('project_',prefix))
project = subset_samples(project, Elevation !=2013)


OTU <- otu_table(project)
TAX <- tax_table(project)
ENV <- sample_data(project)

ENV_dt2 = data.frame(label= rownames(ENV))

ENV_dt2$Fug_shannon=ENV$Shannon
ENV_dt2$Fug_invsimpson=ENV$Invsimpson
ENV_dt2$Fug_richness=ENV$Richness
ENV_dt2$Fug_PD=ENV$PD

ENV_dt = full_join(ENV_dt,ENV_dt2)
write.csv(ENV_dt,file = paste0(folder,'/ENV_dt.csv'))
