library(NST)
library(dplyr)
library(phyloseq)
#################### save dir #######################
folder = '7_stochasticity'
dir.create(folder)
rand.time = 100
nworker = 30
 
load('project_16S_ASV')
#initialize project data
project = get(paste0('project_',prefix)) 
##################### single site check
bad.site = levels(sample_data(project)$Plot)[summary(sample_data(project)$Plot)==1]
project = subset_samples(project,!(Plot %in% bad.site))
#####################
OTU <- otu_table(project) %>% data.frame(.)
ENV <- sample_data(project)
message('data loaded')
treat = ENV[,'Plot'] %>% data.frame(.)
tnstout=tNST(comm=t(OTU), 
             group=treat, 
             dist.method="bray",
             abundance.weighted=TRUE, 
             rand=rand.time,
             output.rand=TRUE,
             nworker=nworker, 
             null.model="PF", 
             between.group=TRUE,
             SES=TRUE, RC=TRUE)
assign(paste0(prefix,'_',note,'_st'),tnstout) 
save(list = c(paste0(prefix,'_',note,'_st'),'prefix','note'),file = paste0(folder,'/',prefix,'_',note,'_st'))

# bootstrapping test for tNST
tnst.bt=NST::nst.boot(nst.result=tnstout, 
                      group=treat,
                      rand=rand.time, 
                      nworker=nworker)
assign(paste0(prefix,'_',note,'_st_boot'),tnst.bt) 
save(list = c(paste0(prefix,'_',note,'_st_boot'),'prefix','note'),file = paste0(folder,'/',prefix,'_',note,'_st_boot'))

########################################
load('project_ITS_ASV')
#initialize project data
project = get(paste0('project_',prefix)) 
note = 'all'
##################### single site check
bad.site = levels(sample_data(project)$Plot)[summary(sample_data(project)$Plot)==1]
project = subset_samples(project,!(Plot %in% bad.site))
#####################
OTU <- otu_table(project) %>% data.frame(.)
ENV <- sample_data(project)
message('data loaded')
treat = ENV[,'Plot'] %>% data.frame(.)
tnstout=tNST(comm=t(OTU), 
             group=treat, 
             dist.method="bray",
             abundance.weighted=TRUE, 
             rand=rand.time,
             output.rand=TRUE,
             nworker=nworker, 
             null.model="PF", 
             between.group=TRUE,
             SES=TRUE, RC=TRUE)
assign(paste0(prefix,'_',note,'_st'),tnstout) 
save(list = c(paste0(prefix,'_',note,'_st'),'prefix','note'),file = paste0(folder,'/',prefix,'_',note,'_st'))

# bootstrapping test for tNST
tnst.bt=NST::nst.boot(nst.result=tnstout, 
                      group=treat,
                      rand=rand.time, 
                      nworker=nworker)
assign(paste0(prefix,'_',note,'_st_boot'),tnst.bt) 
save(list = c(paste0(prefix,'_',note,'_st_boot'),'prefix','note'),file = paste0(folder,'/',prefix,'_',note,'_st_boot'))
