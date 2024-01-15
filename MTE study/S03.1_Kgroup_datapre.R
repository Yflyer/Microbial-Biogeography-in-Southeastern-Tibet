source('0_function.R')
#################### save dir #######################
folder = paste0('S03')
dir.create(folder)
########### input ##############
# set input
Env = read.csv('S01/MTE_ENV.csv',row.names = 1)
load('project_16S_ASV')
project = get(paste0('project_',prefix))
samp =  group_count(otu_table(project),sample_data(project)$Plot)
group_occur = rowSums(samp>0)


Nthreads = 15
Memory.Limitation = 200
null.model = "taxa.labels"
abundance.weighted = FALSE
runs = 99
k_list=c(3:12)

########### data ##############
# Local Scale Data
Phylo_Kgroup_dt = data.frame()
for (k in k_list) {

  project = prune_taxa(taxa_names(project)[group_occur==k], project)
  
  set.seed(1234)
  Min3SampleDepth(project)
  #project = rarefy_even_depth(project,sample.size = min(sample_sums(project)),rngseed = T)
  #project = prune_taxa(taxa_names(project)[sample(1:40000,500)], project)
  
  ######### create big distance matrix
  # Get phylo dist (Big)
  tree = phy_tree(project)
  Save_wd = paste0('PD_',prefix,'_k',k)
  if(!file.exists(file.path(Save_wd,"pd.desc"))) {
    pd.big=pdist.big(tree = tree, 
                     wd=Save_wd, 
                     nworker = 70, 
                     memory.G = Memory.Limitation)
    pd.desc = pd.big$pd.file
    pd.spname = pd.big$tip.label
    pd.wd = pd.big$pd.wd
  }else{
    # if you already calculated the phylogenetic distance matrix in a previous run
    pd.desc = "pd.desc"
    pd.spname = read.csv(file.path(Save_wd,"pd.taxon.name.csv"),row.names = 1,stringsAsFactors = FALSE)[,1]
    pd.wd = Save_wd
    
  }

  
  ####### parallel calculation ############
  ### parallel setting
  if (Nthreads > detectCores()) {  Nthreads = detectCores() }
  sfInit(parallel = T,cpus = Nthreads)
  
  Phylo_dt = data.frame()
  permute.m = shuffleSet(length(pd.spname), runs)
  cat("### Attached a permutation matrix with ",format(object.size(permute.m), units = "auto"),"\n")
  
  # Calculate phylo diversity
    gc()
    
    samp =  otu_table(project) %>% as.matrix(.) %>% t()
    ###  Data container
    # mntd
    NumSample = dim(samp)[1]
    mntd = mpd = numeric(NumSample)
    # mntd null model
    mntd.rand = mpd.rand = numeric(NumSample)
    mntd.rand.m = mpd.rand.m = matrix(nrow = runs,ncol = NumSample,dimnames = list(1:runs,rownames(samp)))
    
    # sub-process loading
    sfExport("samp","pd.wd","pd.desc","pd.spname","abundance.weighted","MT.mntd","MT.mpd","MT.mpd.mntd") # ,"dis"

    mpd.mntd = sfSapply(1:NumSample,function(N) MT.mpd.mntd(N,samp,pd.wd,pd.desc,pd.spname,abundance.weighted))
    mpd= mpd.mntd[1,]
    mntd= mpd.mntd[2,]
    
    cat("### Now start to test permutation null model (",runs,"runs)\n")
    for (i in 1:runs) {
      cat("# Run:",i," --------",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"\n")
      
      pd.spname.rand = pd.spname[permute.m[i,]]
      sfExport("pd.spname.rand")
      mpd.mntd.rand = sfSapply(1:NumSample,function(N) MT.mpd.mntd(N,samp,pd.wd,pd.desc,pd.spname.rand,abundance.weighted))
      mpd.rand= mpd.mntd.rand[1,]
      mntd.rand= mpd.mntd.rand[2,]
      
      mntd.rand.m[i,] = mntd.rand
      mpd.rand.m[i,] = mpd.rand
    }
    
    # generate result
    mntd.rand.mean <- colMeans(mntd.rand.m,na.rm = T)
    mntd.rand.sd <- apply(X = mntd.rand.m, MARGIN = 2, FUN = sd,na.rm = TRUE)
    mntd.z <- (mntd - mntd.rand.mean)/mntd.rand.sd
    
    mpd.rand.mean <- colMeans(mpd.rand.m,na.rm = T)
    mpd.rand.sd <- apply(X = mpd.rand.m, MARGIN = 2, FUN = sd,na.rm = TRUE)
    mpd.z <- (mpd - mpd.rand.mean)/mpd.rand.sd
    
    PD = pd(t(otu_table(project)), phy_tree(project), include.root = TRUE)
    
    # Final result
    Phylo_dt = data.frame(
      Sample_ID = row.names(samp),
      Richness = rowSums(samp>0), 
      Read_count = rowSums(samp),
      PD,
      mntd, 
      mntd.rand.mean, 
      mntd.rand.sd, 
      nri = mntd.z, 
      mpd, 
      mpd.rand.mean, 
      mpd.rand.sd, 
      nti = mpd.z, 
      runs = runs)
  
  sfStop()
  Phylo_dt$K = k
  Phylo_Kgroup_dt = rbind(Phylo_Kgroup_dt,Phylo_dt)
}

# save.result
write.csv(Phylo_Kgroup_dt,file = paste0(folder,'/Phylo_Kgroup_LocalScale_',prefix,'.csv'))

