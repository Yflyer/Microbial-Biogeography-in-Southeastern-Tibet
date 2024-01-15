source('0_function.R')
### pack up fucntion to add index ###
index_packup <- function(project,ENV){
 # alpha diversity index
  ENV$Shannon=diversity(t(otu_table(project)),index='shannon')
  ENV$Simpson=diversity(t(otu_table(project)),index='simpson')
  ENV$Invsimpson=diversity(t(otu_table(project)),index='invsimpson')
  ENV$Richness=specnumber(t(otu_table(project)))
  ENV$Evenness=ENV$Shannon/log(ENV$Richness)
  
  # phylogentic diversity
  PD = pd(t(otu_table(project)), phy_tree(project), include.root = TRUE)
  ENV$PD = PD$PD 

  return(ENV)
}
  

### Meta information  
meta <- read.csv('0_data/env.txt',header = TRUE,row.names = 1,sep='\t')

##############################################
######################################## 16s
prefix = '16S_OTU'
dt <- read.csv('0_data/16S_otu.tsv',header = TRUE,row.names = 1,sep='\t')
taxa <- read.csv('0_data/16S_otu_taxonomy.tsv',header = TRUE,sep='\t',row.names = 1)
tree <- read_tree('0_data/16S_otu_tree.nwk')

################### taxa filter ################
# taxa filter is very important
# 16S bacteria filter
summary(taxa$kingdom)
taxa = taxa[(taxa$kingdom=='Bacteria'),]

### remove species which is only detected once within each plot ###
index = rowSums(dt>0)>1
print(paste0("there are ",sum(!index)," once-detected species"))
dt = dt[index,]

###############################################
##### phylo object
{
  OTU <- otu_table(dt,taxa_are_rows = TRUE) %>% .[,order(colnames(.))]
  TAX <- taxa %>% as.matrix(.) %>% tax_table(.) 
  ENV <- meta  %>% sample_data(.) %>% .[order(rownames(.)),]
}

##################################
########## data process (resample)
set.seed(930827)
colSums(OTU) %>% sort(.)
project = phyloseq(OTU,TAX,ENV,tree) %>% rarefy_even_depth(.,sample.size=40000)

#########
### other indexes
ENV = sample_data(project)
ENV = index_packup(project,ENV)
sample_data(project)=ENV
##########
assign(paste0('project_',prefix),project) 
save(list = c(paste0('project_',prefix),'prefix'),file = paste0('project_',prefix))

### asv
prefix = '16S_ASV'
dt <- read.csv('0_data/16S_asv.tsv',header = TRUE,row.names = 1,sep='\t')
taxa <- read.csv('0_data/16S_asv_taxonomy.tsv',header = TRUE,sep='\t',row.names = 1)
tree <- read_tree('0_data/16S_asv_tree.nwk')

################### taxa filter ################
# taxa filter is very important
# 16S bacteria filter
summary(taxa$kingdom)
taxa = taxa[(taxa$kingdom=='Bacteria'),]

### remove species which is only detected once within each plot ###
index = rowSums(dt>0)>1
print(paste0("there are ",sum(!index)," once-detected species"))
dt = dt[index,]

###############################################
##### phylo object
{
  OTU <- otu_table(dt,taxa_are_rows = TRUE) %>% .[,order(colnames(.))]
  TAX <- taxa %>% as.matrix(.) %>% tax_table(.) 
  ENV <- meta  %>% sample_data(.) %>% .[order(rownames(.)),]
}

##################################
########## data process (resample)
set.seed(930827)
colSums(OTU) %>% sort(.)
project = phyloseq(OTU,TAX,ENV,tree) %>% rarefy_even_depth(.,sample.size=40000)

#########
### other indexes
ENV = sample_data(project)
ENV = index_packup(project,ENV)
sample_data(project)=ENV
##########
assign(paste0('project_',prefix),project) 
save(list = c(paste0('project_',prefix),'prefix'),file = paste0('project_',prefix)) 

################################################
#############################   ITS_ASV ########
prefix = 'ITS_ASV'
dt <- read.csv('0_data/ITS_asv.tsv',header = TRUE,row.names = 1,sep='\t') 
taxa <- read.csv('0_data/ITS_asv_taxonomy.tsv',header = TRUE,sep='\t',row.names = 1)
tree <- read_tree('0_data/ITS_asv_tree.nwk')

################### taxa filter ################
# taxa filter is very important
# 16S bacteria filter
summary(as.factor(taxa$kingdom))
taxa = taxa[(taxa$kingdom=='Fungi'),]

### remove species which is only detected once within each plot ###
index = rowSums(dt>0)>1
print(paste0("there are ",sum(!index)," once-detected species"))
dt = dt[index,]

###############################################
##### phylo object
{
  OTU <- otu_table(dt,taxa_are_rows = TRUE) %>% .[,order(colnames(.))]
  TAX <- taxa %>% as.matrix(.) %>% tax_table(.) 
  ENV <- meta  %>% sample_data(.) %>% .[order(rownames(.)),]
}

##################################
########## data process (resample)
set.seed(930827)
colSums(OTU) %>% sort(.)
project = phyloseq(OTU,TAX,ENV,tree) %>% rarefy_even_depth(.,sample.size=40000)

#########
### other indexes
ENV = sample_data(project)
ENV = index_packup(project,ENV)
sample_data(project)=ENV
##########
assign(paste0('project_',prefix),project) 
save(list = c(paste0('project_',prefix),'prefix'),file = paste0('project_',prefix)) 

######################################## 
prefix = 'ITS_OTU'
dt <- read.csv('0_data/ITS_otu.tsv',header = TRUE,row.names = 1,sep='\t') 
taxa <- read.csv('0_data/ITS_otu_taxonomy.tsv',header = TRUE,sep='\t',row.names = 1)
tree <- read_tree('0_data/ITS_otu_tree.nwk')

################### taxa filter ################
# taxa filter is very important
# 16S bacteria filter
summary(as.factor(taxa$kingdom))
taxa = taxa[(taxa$kingdom=='Fungi'),]

### remove species which is only detected once within each plot ###
index = rowSums(dt>0)>1
print(paste0("there are ",sum(!index)," once-detected species"))
dt = dt[index,]

###############################################
##### phylo object
{
  OTU <- otu_table(dt,taxa_are_rows = TRUE) %>% .[,order(colnames(.))]
  TAX <- taxa %>% as.matrix(.) %>% tax_table(.) 
  ENV <- meta  %>% sample_data(.) %>% .[order(rownames(.)),]
}

##################################
########## data process (resample)
set.seed(930827)
colSums(OTU) %>% sort(.)
project = phyloseq(OTU,TAX,ENV,tree) %>% rarefy_even_depth(.,sample.size=40000)

#########
### other indexes
ENV = sample_data(project)
ENV = index_packup(project,ENV)
sample_data(project)=ENV
##########
assign(paste0('project_',prefix),project) 
save(list = c(paste0('project_',prefix),'prefix'),file = paste0('project_',prefix)) 

