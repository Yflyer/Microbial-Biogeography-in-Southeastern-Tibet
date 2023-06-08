source('0_function.R')

#################### save dir #######################
folder = paste0('3_disper')
dir.create(folder)

##################################################
### Select the project that you require
load('project_16S_OTU')
load('project_16S_ASV')
load('project_ITS_OTU')
load('project_ITS_ASV')

##################################################
#initialize project data
project = get(paste0('project_',prefix)) 
OTU <- otu_table(project) %>% data.frame()
ENV <- sample_data(project) %>% data.frame()
Dt_nmds = vegdist(t(OTU)) %>% metaMDS(.,k=3)

# Examine NMDS
Dt_nmds$stress # <0.2
Dt_nmds_site <- data.frame(Dt_nmds$points)

# Check names
rownames(ENV) == rownames(Dt_nmds_site)

plot_dt=cbind(ENV,Dt_nmds_site)

# get points of convex hull
convex.dt = data.frame()
ENV$Plot=as.factor(ENV$Plot)
for (i in 1:nlevels(ENV$Plot)) {
  hull.point = Dt_nmds_site[ENV$Plot==levels(ENV$Plot)[i],1:2] %>% chull(.) 
  conv.hull.subdt = plot_dt[ENV$Plot==levels(ENV$Plot)[i],] %>% .[hull.point,]
  convex.dt = rbind(convex.dt,conv.hull.subdt)
}

##################################################
# Plot
PTs = theme(text = element_text(size=16),
            panel.background=element_blank(),
            legend.position='none',
            panel.border = element_rect(colour = "grey", fill=NA, linewidth=1))
ggplot(data=plot_dt)+
  geom_point(aes(x=MDS1,y=MDS2,size=normalization(Rs),color=Elevation),alpha=0.4)+
  geom_polygon (data=convex.dt,aes(x=MDS1,y=MDS2,group=Plot,color=Elevation),size=0.5,fill=NA)+
  PTs+
  scale_color_viridis()+
  ylim(-0.5,0.5)+
  xlim(-0.5,0.5)
ggsave(filename = paste0(folder,'/',prefix,'_dispersal_hull.pdf'),dpi=900,width=11,height=14,units='cm')



