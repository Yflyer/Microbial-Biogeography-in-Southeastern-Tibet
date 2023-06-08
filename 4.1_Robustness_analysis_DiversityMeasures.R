source('0_function.R')
########################
folder = paste0('4_Robustness_analysis')
dir.create(folder)

########################################
### Define function: Diversity measures along elevation by different phylum
get_measures_data <-function(project,prefix,measures_dt,Elevation,taxa_list,phylum_list){
  plot_data = data.frame()
  for (i in 0:length(phylum_list)) {
    # set overall information
    if (i==0) {    
      sub_dt = otu_table(project)
    } else {
      sub_dt = otu_table(project)[taxa_list==phylum_list[i],] 
    }
    
    count_dt = colSums(sub_dt)
    # Ensure phylum could be detected in all samples
    if (all(count_dt>0)) {
      measures_dt$Abundance = colSums(sub_dt)
      ### alpha
      measures_dt$Richness = specnumber(t(sub_dt))
      measures_dt$Shannon = diversity(t(sub_dt),index='shannon')
      measures_dt$Evenness = measures_dt$Shannon/measures_dt$Richness
      measures_dt$Inv.Simpson  = diversity(t(sub_dt),index='invsimpson')
      
      ### beta (disper)
      # bray 
      Beta = vegdist(t(sub_dt), method="bray", binary=FALSE,na.rm = T)
      measures_dt$Bray.Curtis = betadisper(Beta,Elevation)$distances
      
      # Sor
      Beta = vegdist(t(sub_dt), method="bray", binary=T,na.rm = T)
      measures_dt$Sorensen = betadisper(Beta,Elevation)$distances
      
      # horn
      Beta = vegdist(t(sub_dt), method="horn", binary=T,na.rm = T)
      measures_dt$Horn = betadisper(Beta,Elevation)$distances
      
      # jaccard
      Beta = vegdist(t(sub_dt), method="jaccard", binary=T,na.rm = T)
      measures_dt$Jaccard = betadisper(Beta,Elevation)$distances
      
      ### other label
      # set overall information
      if (i==0) {        
        measures_dt$Taxa_name = paste0('Overall community (',prefix,')')
        measures_dt$Level = 'Overall'
      } else {
        measures_dt$Taxa_name = phylum_list[i]
        measures_dt$Level = 'Phylum'
      }
      
      measures_dt$Data = prefix
      plot_data = rbind(plot_data,measures_dt)
    }
  }
  return(plot_data)
}

############### Env data pre
ENV <- read.csv('1_Data_Merge/Merge_dt.csv',row.names = 1,stringsAsFactors = T)
Elevation = ENV$Elevation %>% as.factor()
measures_dt = select(ENV,Sample_ID,Elevation,Ecosystem_type,Plot,Rs,SOC,TN,Temp.)
######################## 16S
load('project_16S_ASV')
project=get(paste0('project_',prefix))
### generate candidate of bacteria
TAX <- tax_table(project)
Level = "phylum"

### taxa list for level
taxa_list = TAX[,Level]
taxa_list[taxa_list == ""] = "unknownbacteria"
OTU <- otu_table(project) %>% group_sum(.,taxa_list,margin = 2)

### select number of Taxa by mean abundance (remove zero) > 0.1% at phylum level
RA_phylum = rowMeans(OTU,na.rm = T)  %>% sort(.,decreasing = T)
round(RA_phylum/40000,3)
phylum_list = names(RA_phylum)[RA_phylum> (40000 * 0.001)]

### Get measures data of bacteria
bac_data = get_measures_data(project,prefix,measures_dt,Elevation,taxa_list,phylum_list)


######################## ITS
load('project_ITS_ASV')
project=get(paste0('project_',prefix))
### generate candidate of fungi
TAX <- tax_table(project)
Level = "phylum"

### taxa list for level
taxa_list = TAX[,Level]
taxa_list[taxa_list == ""] = "unknownfungi"
OTU <- otu_table(project) %>% group_sum(.,taxa_list,margin = 2)

### select number of Taxa by mean abundance (remove zero) > 0.1% at phylum level
RA_phylum = rowMeans(OTU,na.rm = T)  %>% sort(.,decreasing = T)
round(RA_phylum/40000,3)

phylum_list = names(RA_phylum)[RA_phylum> (40000 * 0.001)]

### Get measures data of fungi
fug_data = get_measures_data(project,prefix,measures_dt,Elevation,taxa_list,phylum_list)

### merge
final_data = rbind(bac_data,fug_data)

write.csv(final_data,paste0(folder,'/Taxa_diversity_data.csv'))
