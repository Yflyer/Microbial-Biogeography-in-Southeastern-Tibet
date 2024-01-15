source('0_function.R')
library(ggplot2)
library(ggh4x)
# library(spaa)
#################### save dir #######################
folder = paste0('SX')
dir.create(folder)
########### data ##############
# load data
Env = read.csv('S02/MTE_ENV.csv',row.names = 1,stringsAsFactors = T)
Process_bin = read.csv('S01/icamp_16S_ASV/16S_ASV.BinContributeToProcess_EachGroup.csv')
Taxa_Bin = read.csv('S01/icamp_16S_ASV/16S_ASV.Taxon_Bin.csv')
Rep_Bin = read.csv('S01/icamp_16S_ASV/16S_ASV.Bin_TopTaxon.csv')

load('project_16S_ASV')
project = get(paste0('project_',prefix))
project = prune_taxa(Rep_Bin$TopTaxonID,project)
OTU = otu_table(project) %>% data.frame() #%>% .[Taxa_Bin$ID,]

# K data
samp =  group_count(OTU,sample_data(project)$Plot)
K = rowSums(samp>0)
Occurrence = rowSums(OTU>0)

# Niche breadth
Niche.breadth = niche.width(t(OTU), method = c( "levins")) %>% as.numeric()

# stochastic data
Sto_dt = Rep_Bin %>% mutate(K, Occurrence,Niche.breadth) %>%
  group_by(Bin) %>%
  summarise(K = mean(K),
            Occurrence = mean(Occurrence),
            Richness = n(),
            Niche.breadth = mean(Niche.breadth))

dt = data.frame(Bin = Rep_Bin$Bin,
                label = Rep_Bin$TopTaxonID,
                Relative.abundance = Rep_Bin$TopTaxonRAinBin,
                Phylum = Rep_Bin$TopTaxon.phylum,
                Class = Rep_Bin$TopTaxon.class)
Tax_dt = Sto_dt %>% 
  left_join(dt) %>% 
  mutate(Phylum=gsub('_','',Phylum),
         Class=gsub('_','',Class)) %>%
  filter(Class !='')


Order.taxa = Tax_dt$Phylum[order(Tax_dt$Relative.abundance,decreasing = T)] %>% unique(.) %>% .[1:10]
Tax_dt$Phylum[!Tax_dt$Phylum %in% Order.taxa]='Others'
Tax_dt$Phylum = factor(Tax_dt$Phylum,levels=c(Order.taxa,'Others'))
rownames(Tax_dt) = Tax_dt$label

# Bin process
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

#
library(ggtree)
library(ggtreeExtra)
library(treeio)
tree <- prune_taxa(Tax_dt$label,project) %>% phy_tree() %>% as.treedata()
#tree <- full_join(as_tibble(tree), Tax_dt, by = 'label') %>% as.treedata()
data=fortify(tree)

# basic tree plot
p <- ggtree(tree,size=0.4,color='grey50')

Tax_dt = Tax_dt %>% relocate(label) %>% as_tibble()
p2 <- p +
  geom_fruit(
    data=Tax_dt,
    geom=geom_point,
    mapping=aes(y=label, color=Phylum),
    size=1,
    #shape=21,
    #fill='white',# , size=Genome_size
    position="identity"
  )+
  scale_size_continuous(
    range=c(0.05, 1.5), # the range of size.
    guide=guide_legend(
      keywidth=1, 
      keyheight=1)
  )
p2

# Merge
# Bin process
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

process_dt = Tax_dt %>%
  left_join(Between_dt) %>%
  drop_na()

p3 <- p2 + 
  geom_fruit(
    data=process_dt,
    geom=geom_tile,
    fill = 'steelblue',
    mapping=aes(y=label, x=Process, alpha=ST),
    pwidth=0.2,
    axis.params=list(
      axis="x", # add axis text of the layer.
      text.angle=-45,
      text.size=3,# the text angle of x-axis.
      hjust=1  # adjust the horizontal position of text of axis.
    )
  )

### Mantel r data
Bin_list = Tax_dt$Bin %>% unique()
library(vegan)
load('project_16S_ASV')
project = get(paste0('project_',prefix))
bin_mantel = function(Bin_ID,Taxa_Bin,project){
  require(dplyr)
  require(phyloseq)
  require(vegan)
  sub_taxa = Taxa_Bin$ID[Taxa_Bin$Bin == Bin_ID]
  sub_project = prune_taxa(sub_taxa,project)
  sub_OTU = otu_table(sub_project) 
  # remove empty sample
  sub_OTU = sub_OTU[,colSums(sub_OTU>0)>0]
  Sample = colnames(sub_OTU)
  Sub_Env = Env %>% .[Sample,]
  
  # generate dist
  Bray = sub_OTU %>% t() %>% vegdist() # binary = T
  Temp = Sub_Env$Temp. ; names(Temp) = rownames(Sub_Env)
  
  # Richness
  sub_Richness = specnumber(t(sub_OTU))
  Cor.Result = cor.test(sub_Richness,Temp)
  
  Other = Sub_Env %>% data.frame()%>% select(SOC,TN,Water_content,pH) %>% scale() %>% dist()
  Temp = Temp %>% scale() %>% dist()
  Man.Result = mantel.partial(Bray,Temp,Other, method="spearman")
  #Man.c.Result = mantel(Bray,Temp, method="spearman")
  
  # get result
  v = c(Bin = Bin_ID,
        Richness = length(sub_taxa),
        Cor = as.numeric(Cor.Result$estimate),
        Cor.P = Cor.Result$p.value,
        Mantelr = Man.Result$statistic,
        Sign = Man.Result$signif)
  return(v)
}
library(snowfall)
sfInit(parallel=TRUE, cpus=12)
sfExport( "Env")
Result = sfSapply(Bin_list,bin_mantel,Taxa_Bin,project)
sfStop()
mantel_dt = t(Result) %>% data.frame() %>%
  left_join(Tax_dt[,c('Bin','label','K','Phylum')])
mantel_dt$Cor=as.numeric(mantel_dt$Cor)
mantel_dt$Cor.P=as.numeric(mantel_dt$Cor.P)
mantel_dt$Mantelr=as.numeric(mantel_dt$Mantelr)
mantel_dt$Sign=as.numeric(mantel_dt$Sign)
mantel_dt$P='P <= 0.001'
mantel_dt$P[mantel_dt$Sign>0.010]='P < 0.050'
mantel_dt$P[mantel_dt$Sign>0.05]='P > 0.050'
mantel_dt$K_group='K<8'
mantel_dt$K_group[mantel_dt$K>=9]='K>=8'
p4 <- p3 +  
  geom_fruit(
    data=mantel_dt,
    geom=geom_bar,
    mapping=aes(y=label, x=Mantelr, fill=K),  #The 'Cazy_num' of 'dat1' will be mapped to x
    pwidth=0.2,
    stat="identity",
    orientation="y", # the orientation of axis.
    axis.params=list(
      axis="x", # add axis text of the layer.
      text.angle=-45, # the text size of axis.
      hjust=0  # adjust the horizontal position of text of axis.
    ),
    grid.params=list() # add the grid line of the external bar plot.
  )
p4
ggsave(filename = paste0(folder,'/','Bin_tree.pdf'),dpi=600,width=16,height=22,units='cm')

ggplot(mantel_dt)+geom_point(aes(K,Cor))
ggplot(mantel_dt)+geom_point(aes(K,Mantelr))

