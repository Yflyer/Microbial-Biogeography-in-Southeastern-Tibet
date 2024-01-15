source('0_function.R')
library(iCAMP)
library(NST)
library(ieggr)
library(tidyr)
cat("##############################\n")


#####################################################
load('project_16S_ASV')
#load('project_ITS_ASV')
#initialize project data
project = get(paste0('project_',prefix)) 
cat("### start project", prefix," \n")
#################### save dir #######################
folder = paste0('7_icamp_',prefix)
dir.create(folder)

# phyloseq transformation for test
#project = subset_samples(project,Plot %in% c('L1','L2'))
OTU <- otu_table(project) %>% data.frame(.) 
#OTU = OTU[sample(1:nrow(OTU),5000),]

# at least 6 times of ASV for correlation
OTU = filter(OTU, rowSums(OTU>0) > 6) %>% drop_na() 
cat("### ", nrow(OTU)," OTUs are remained after fitler\n")

OTU = otu_table(OTU,taxa_are_rows = TRUE)
otu_table(project) = OTU

# extract matched components of phyloseq
OTU <- otu_table(project) %>% data.frame(.) 
TAX <- tax_table(project) %>% data.frame(.)
ENV <- sample_data(project) %>% data.frame(.)
Tre <- phy_tree(project)
treat=ENV %>% mutate(Elevation=as.factor(Elevation)) %>% select(Elevation)
Var = ENV %>% select(Rs) # used to filter ENV to construct niche 
save(list=c('project','prefix','OTU','TAX','ENV','Tre','treat','Var'),file =paste0(folder,'/icamp_preload_',prefix,'.rda'))

###############################################
##############################