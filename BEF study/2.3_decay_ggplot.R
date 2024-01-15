source('0_function.R')
########################
folder = paste0('2_decay')
dir.create(folder)
######################## load data
ENV <- read.csv('1_Data_Merge/Merge_dt.csv',row.names = 1,stringsAsFactors = T)
decay_data = read.csv(file = paste0(folder,'/decay_data.csv'))

decay_data$Bac.Similarity = 1-decay_data$Bac_Bray
decay_data$Fug.Similarity = 1-decay_data$Fug_Bray
decay_data$Log.distance = log(decay_data$Distance+1,10)

### overall scale
PTs <- theme(text = element_text(size=12),
             panel.background=element_blank(),
             legend.position = "none",
             panel.border = element_rect(colour = "grey", fill=NA, linewidth=1))
################################################################
### bac
plot.data = decay_data

SDecay <- ggplot(plot.data,aes(x=Log.distance,y=Bac.Similarity))+
  geom_point(alpha = 0.2,stroke = 0.1,size=1,shape=16,color='grey20')
Lm_overall = geom_smooth(method='lm',data = plot.data,aes(x=Log.distance,y=Bac.Similarity),color='gray20',se=F,linewidth=1)
SDecay = SDecay+Lm_overall +PTs
SDecay

### overall
model.dt = plot.data
### distance slope
decay <- lm(Bac.Similarity ~ Log.distance, data=model.dt)
lm.slope = summary(decay)$coefficients %>% .['Log.distance','Estimate'] %>% specify_decimal(.,3)
lm.P = lmp(decay) %>% labelpstar()
Overall_Slope_text = annotate("text", y=0.85,x=1,label= paste('Slope(Dis.): ',lm.slope,lm.P),)

SDecay+Overall_Slope_text
ggsave(filename = paste0(folder,'/decay_bac.bray.pdf'),dpi=900,width=8,height=7,units='cm')

################################################################
### fungal
plot.data = decay_data
SDecay <- ggplot(plot.data,aes(x=Log.distance,y=Fug.Similarity))+
  geom_point(alpha = 0.2,stroke = 0.1,size=1,shape=16,color='grey20') +PTs

Lm_overall = geom_smooth(method='lm',data = plot.data,aes(x=Log.distance,y=Fug.Similarity),color='gray20',se=F,linewidth=1)
SDecay = SDecay+Lm_overall

### overall
model.dt = plot.data
### distance slope
decay <- lm(Fug.Similarity ~ Log.distance, data=model.dt)
lm.slope = summary(decay)$coefficients %>% .['Log.distance','Estimate'] %>% specify_decimal(.,3)
lm.P = lmp(decay) %>% labelpstar()
Overall_Slope_text = annotate("text", y=0.85,x=1,label= paste('Slope(Dis.): ',lm.slope,lm.P))

SDecay+Overall_Slope_text
ggsave(filename = paste0(folder,'/decay_Fug.bray.pdf'),dpi=900,width=8,height=7,units='cm')

