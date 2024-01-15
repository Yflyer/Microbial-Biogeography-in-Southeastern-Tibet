source('0_function.R')
########################
folder = paste0('4_Robustness_analysis')
dir.create(folder)
########################
### load plot data
plot_data = read.csv(paste0(folder,'/Taxa_diversity_data.csv'),row.names = 1)

### prepare x axis level
bac_index = plot_data$Taxa_name[plot_data$Data=='16S_ASV'] %>% unique() 
fug_index = plot_data$Taxa_name[plot_data$Data=='ITS_ASV'] %>% unique() 
tax_index = c(as.character(bac_index),as.character(fug_index)) 
plot_data$Taxa_name = factor(plot_data$Taxa_name,levels = tax_index  )


### bubble data
colnames(plot_data)
indice = colnames(plot_data)[9:21]
indice = c("Abundance","Richness","Evenness","Shannon","Inv.Simpson","Bray.Curtis","Jaccard","Sorensen","Horn")
bubble_dt = data.frame()

for (i in 1:length(indice)) {
  lm.P =lm.Slope =c()
  for (j in 1:length(tax_index)) {
    sub_dt =  plot_data[plot_data$Taxa_name==tax_index[j],]
    lm_result = cor.test(sub_dt$Rs,sub_dt[[indice[i]]],method = 'spearman')
    lm.P = c(lm.P,lm_result$p.value)
    lm.Slope = c(lm.Slope,lm_result$estimate)
  }
  sub_dt = plot_data %>% group_by(Taxa_name) %>% summarise(Level=unique(Level))
  sub_dt$lm.P = lm.P
  sub_dt$lm.Slope = lm.Slope
  sub_dt$indice = indice[i]
  bubble_dt=rbind(bubble_dt,sub_dt)
}

# stat 
bubble_dt$lm.P = p.adjust(bubble_dt$lm.P, method = 'BH')
bubble_dt$lm.P[is.na(bubble_dt$lm.P)] = 1
bubble_dt$Star = sapply(bubble_dt$lm.P,function(x)labelpstar(x))
bubble_dt$Correlation = 'Uncorrelated'
bubble_dt$Correlation[ bubble_dt$lm.P<0.05] = 'Correlated'
bubble_dt$indice = factor(bubble_dt$indice,levels = indice)

write.csv(bubble_dt,paste0(folder,'/correlation_data.csv'))

######################################################################
# rephrase slope data to output as a table
slope_dt = bubble_dt %>% 
  select(Taxa_name,indice,lm.Slope) %>%
  mutate(lm.Slope = round(lm.Slope,3)) %>%
  spread(indice,lm.Slope)
slope_dt[is.na(slope_dt)]=0

p_dt = bubble_dt %>% 
  select(Taxa_name,indice,lm.P) %>%
  spread(indice,lm.P)

slope_dt[p_dt>=0.05] = NA

output=data.frame()
for (i in 1:nrow(slope_dt)) {
  output_row = paste0(slope_dt[i,-1],sapply(p_dt[i,-1], labelpstar))
  output=rbind(output,output_row)
}
colnames(output) = colnames(slope_dt)[-1]
rownames(output) = slope_dt$Taxa_name

write.csv(output,paste0(folder,'/slope_dt.csv'))

####################################
### bubble plot
PTs <- theme(text = element_text(size=10),
             axis.text.x = element_text(hjust = 0.5,vjust = 0.5,angle = 90),
             panel.background=element_blank(),
             legend.position = "bottom",
             panel.border = element_rect(colour = "black", fill=NA, size=1))
sfg = scale_fill_gradient2(low="#6E9BC5", mid="white", high="#F091A0", #colors in the scale
                           midpoint=0,    
                           breaks=seq(-0.5,0.5,0.25), #breaks in the scale bar
                           limits=c(-0.6, 0.6) #same limits for plots
                           )

ggplot(bubble_dt)+geom_tile(aes(x=Taxa_name,y=indice,fill=lm.Slope,alpha=Correlation))+scale_alpha_discrete(range=c(1, 0))+sfg+PTs
ggsave(filename = paste0(folder,'/Indice_bubble_plot.pdf'),dpi=900,width=17,height=11,units='cm')

