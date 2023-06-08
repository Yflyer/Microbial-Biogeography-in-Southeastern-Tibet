source('0_function.R')
########################
folder = paste0('2_decay')
dir.create(folder)
##################################################
###### Get environment distance matrix 
ENV <- read.csv('1_Data_Merge/Merge_dt.csv',row.names = 1,stringsAsFactors = T)

### create within-group distance data
WG.distance = dist(ENV[,c('v_axis','h_axis')]) %>% col3m(.,dist_name = 'WG.distance')
Env_data = WG.distance

### create environmental factors
colnames(ENV)
Spatialfactor = c('Temp.','pH') 
EcoFun = c('Rs')
Nutrient = c('Water_content','SOC','TN','NO3','NH4')
Var=c(Spatialfactor,EcoFun,Nutrient) 

### Get distance matrix by Euclidean distance
for (i in 1:length(Var)) {
  factor_dt = dist(ENV[,Var[i], drop = F],upper=T) %>% col3m(.,dist_name = Var[i]) 
  Env_data = inner_join(Env_data,factor_dt)
}

### add other label as factors
DiscreFactor = c('Ecosystem_type','Plot','Elevation')
for (i in 1:length(DiscreFactor)) {
  row = sapply(Env_data$row,function(x) ENV[[DiscreFactor[i]]][rownames(ENV)==x]) %>% as.character()
  col = sapply(Env_data$col,function(x) ENV[[DiscreFactor[i]]][rownames(ENV)==x]) %>% as.character()
  Env_data[[paste0(DiscreFactor[i],'_row')]]=row
  Env_data[[paste0(DiscreFactor[i],'_col')]]=col
}

Env_data$Scale='Regional'
Env_data$Scale[Env_data$Plot_col==Env_data$Plot_row]='Within-group'
Env_data$WG.distance[Env_data$Scale!='Within-group']=0
Distance = geo_dist(ENV,Lon = 'Lon.',Lat = 'Lat.') %>% col3m(.,dist_name = 'Distance')
Env_data=inner_join(Env_data,Distance,factor_dt,by=c('row','col'))
Env_data$Distance = Env_data$Distance + Env_data$WG.distance * 0.001

###### Get beta distance matrix by multiple measures
### load project 16S
load('project_16S_ASV')
load('project_ITS_ASV')
### Bray.Curtis distance
Bac_Bray = otu_table(project_16S_ASV) %>% t(.) %>% vegdist(., method="bray", binary=FALSE) %>% col3m(.,dist_name ='Bac_Bray')
Bac_Sorensen = otu_table(project_16S_ASV) %>% t(.) %>% vegdist(., method="bray", binary=TRUE) %>% col3m(.,dist_name ='Bac_Sorensen')
Fug_Bray = otu_table(project_ITS_ASV) %>% t(.) %>% vegdist(., method="bray", binary=FALSE) %>% col3m(.,dist_name ='Fug_Bray')
Fug_Sorensen = otu_table(project_ITS_ASV) %>% t(.) %>% vegdist(., method="bray", binary=TRUE) %>% col3m(.,dist_name ='Fug_Sorensen')

### Phylogentic distance
# Get binary tree before unifrac
library(ape)
tree <- phy_tree(project_16S_ASV)
is.binary(tree)
phy_tree(project_16S_ASV) = multi2di(tree)

tree <- phy_tree(project_ITS_ASV)
is.binary(tree)
phy_tree(project_16S_ASV) = multi2di(tree)

# register parallel back-end
library(doParallel)
cl = makeCluster(detectCores())
registerDoParallel(cl)

# Unifrac distance
Bac_Unifrac = UniFrac(project_16S_ASV, weighted=T, normalized=TRUE, parallel=T, fast=TRUE) %>% col3m(.,dist_name = 'Unifrac')
Bac_Uw.Unifrac = UniFrac(project_16S_ASV, weighted=F, normalized=TRUE, parallel=T, fast=TRUE) %>% col3m(.,dist_name = 'Uw.Unifrac')
Fug_Unifrac = UniFrac(project_16S_ASV, weighted=T, normalized=TRUE, parallel=T, fast=TRUE) %>% col3m(.,dist_name = 'Unifrac')
Fug_Uw.Unifrac = UniFrac(project_16S_ASV, weighted=F, normalized=TRUE, parallel=T, fast=TRUE) %>% col3m(.,dist_name = 'Uw.Unifrac')

stopCluster(cl)

##################################################
### merge distance
decay_data = Reduce(function(x, y) inner_join(x, y,by=c('row','col')),
                    list(
                      Env_data,
                      Bac_Unifrac,Bac_Uw.Unifrac,
                      Bac_Bray,Bac_Sorensen,
                      Fug_Unifrac,Fug_Uw.Unifrac,
                      Fug_Bray,Fug_Sorensen))

### assign data
write.csv(decay_data,paste0(folder,'/decay_data.csv'))

