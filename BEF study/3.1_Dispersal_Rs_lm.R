source('0_function.R')
#################### save dir #######################
folder = paste0('3_Disper_analysis')
dir.create(folder)
####################
plot_data = read.csv('1_Data_Merge/Merge_dt.csv',row.names = 1,stringsAsFactors = T)

cor.test(plot_data$Bac_disper,plot_data$Rs,method = c("pearson"))
cor.test(plot_data$Fug_disper,plot_data$Rs,method = c("pearson"))

lm(Rs ~ Bac_disper,data = data.frame(plot_data)) %>% summary()
lm(Rs ~ Fug_disper,data = data.frame(plot_data)) %>% summary()

PTs <- theme(text = element_text(size=20),
             panel.background=element_blank(),
             legend.position = "none",
             panel.border = element_rect(colour = "black", fill=NA, linewidth=1))

########### Correlation between dispersal and Rs
# bac
ggplot(plot_data,aes(x=Bac_disper,y=Rs))+geom_point(aes(color=Elevation))+PTs+geom_smooth(method='lm',color='black',se=F,linewidth=1.5)+scale_color_viridis()
ggsave(filename = paste0(folder,'/Rs_Bac.Dispersal.pdf'),dpi=900,width=9,height=7,units='cm')

# fug
ggplot(plot_data,aes(x=Fug_disper,y=Rs))+geom_point(aes(color=Elevation))+PTs+geom_smooth(method='lm',color='black',se=F,linewidth=1.5)+scale_color_viridis()
ggsave(filename = paste0(folder,'/Rs_Fug.Dispersal.pdf'),dpi=900,width=9,height=7,units='cm')
colnames(plot_data)

### Another pattern by pooling data within each elevation
pool_data = plot_data %>%
  group_by(Elevation) %>%
  summarise(Bac_disper.avg = mean(Bac_disper),
            Bac_disper.sd = sd(Bac_disper),
            Fug_disper.avg = mean(Fug_disper),
            Fug_disper.sd = sd(Fug_disper),
            Rs.avg = mean(Rs),
            Rs.sd = sd(Rs))

# Average correlation
ggplot(pool_data,aes(x=Fug_disper.avg,y=Rs.avg,color=Elevation))+
  geom_errorbar(aes(ymax=Rs.avg+Rs.sd,ymin=Rs.avg-Rs.sd))+
  geom_errorbar(aes(xmax=Fug_disper.avg+Fug_disper.sd,xmin=Fug_disper.avg-Fug_disper.sd))+
  geom_point(size=3)+
  PTs+geom_smooth(method='lm',color='black',se=F,size=1)+scale_color_viridis()
lm(Rs.avg~Fug_disper.avg,pool_data) %>% summary()

ggplot(pool_data,aes(x=Bac_disper.avg,y=Rs.avg,color=Elevation))+
  geom_errorbar(aes(ymax=Rs.avg+Rs.sd,ymin=Rs.avg-Rs.sd))+
  geom_errorbar(aes(xmax=Bac_disper.avg+Bac_disper.sd,xmin=Bac_disper.avg-Bac_disper.sd))+
  geom_point(size=3)+
  PTs+geom_smooth(method='lm',color='black',se=F,size=1)+scale_color_viridis()
lm(Rs.avg~Bac_disper.avg,pool_data) %>% summary()

# deviation correlation
ggplot(pool_data,aes(x=Fug_disper.sd,y=Rs.sd))+
  geom_point(aes(color=Elevation))+
  PTs+geom_smooth(method='lm',color='black',se=F,size=1.5)+scale_color_viridis()
lm(Rs.sd~Fug_disper.sd,pool_data) %>% summary()

ggplot(pool_data,aes(x=Bac_disper.sd,y=Rs.sd))+
  geom_point(aes(color=Elevation))+
  PTs+geom_smooth(method='lm',color='black',se=F,size=1.5)+scale_color_viridis()
lm(Rs.sd~Bac_disper.sd,pool_data) %>% summary()
