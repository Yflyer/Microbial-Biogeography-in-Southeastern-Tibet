source('0_function.R')
#################### save dir #######################
folder = paste0('S02')
dir.create(folder)
######################
library("doParallel")
library("foreach")

cl = makeCluster(detectCores() ) 
registerDoParallel(cl) 

load('project_16S_ASV')
project = get(paste0('project_',prefix)) 
OTU <- otu_table(project)
Soil.Bac.phylo = UniFrac(project,weighted = F,parallel = T)
Soil.Bac.manhattan = vegdist(t(OTU), method="manhattan", binary=FALSE) 
Soil.Bac.bray = vegdist(t(OTU), method="bray", binary=FALSE)
Soil.Bac.sorensen = vegdist(t(OTU), method="bray", binary=TRUE)
Soil.Bac.Jaccard = vegdist(t(OTU), method="Jaccard", binary=TRUE)

load('project_ITS_ASV')
project = get(paste0('project_',prefix)) 
OTU <- otu_table(project)
Soil.Fug.phylo = UniFrac(project,weighted = F,parallel = T)
Soil.Fug.manhattan = vegdist(t(OTU), method="manhattan", binary=FALSE) 
Soil.Fug.bray = vegdist(t(OTU), method="bray", binary=FALSE)
Soil.Fug.sorensen = vegdist(t(OTU), method="bray", binary=TRUE)
Soil.Fug.Jaccard = vegdist(t(OTU), method="Jaccard", binary=TRUE)

stopCluster()
class(Stream.Fug.sorensen)
Dist.list = list(Soil.Bac.phylo,Soil.Bac.manhattan,Soil.Bac.bray,Soil.Bac.sorensen,Soil.Bac.Jaccard,
                 Soil.Fug.phylo,Soil.Fug.manhattan,Soil.Fug.bray,Soil.Fug.sorensen,Soil.Fug.Jaccard)
names(Dist.list) = c('Soil.Bac.phylo','Soil.Bac.manhattan','Soil.Bac.bray','Soil.Bac.sorensen','Soil.Bac.Jaccard',
                        'Soil.Fug.phylo','Soil.Fug.manhattan','Soil.Fug.bray','Soil.Fug.sorensen','Soil.Fug.Jaccard')
saveRDS(Dist.list,paste0(folder,'/Dist.list'))
