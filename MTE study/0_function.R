############## Yufei's toolkit ###################
# Author: Yufei Zeng
# Email: yfzeng0827@hotmail.com
# web: https://github.com/Yflyer

### load required package
# data tool
library(dplyr)
library(tidyr)
library(scales)
library(stringr)
library(usedist)
library(phyloseq)
library(geosphere)
library(vroom)

# plot
library(ggplot2)
library(ggh4x)
library(ggsci)
#library(ggtree)
library(ggthemes)
library(RColorBrewer)
library(viridis)

# ecology
library(vegan)
library(picante)
library(NST)
library(iCAMP)

### Subjects explanation
# table: rownames are features (OTUs, ASVs, or other  features); colnames are samples
# map: rownames are samples; colnames are features of sample as metadata
# group: a matched treatment information, grouping information or other categorized information could be used for table along samples. Generally extracted from map.

### group operation tools
### clean NA and cut data by group

############ function used #############################
group_count<-function(table,group){
  ### this function is used to count the occurrence frequency of each row in a table in a specific grouping condition
  ### such as: group_count(OTU,treatment)
  ### then we will get each OTU occurrence assigned on treatment levels on rows
  group= group %>% as.factor(.) %>% droplevels(.) # make factor
  count_table = matrix(ncol= nlevels(group),nrow = nrow(table),dimnames = list(rownames(table),levels(group)))
  for (i in 1:ncol(count_table)) {
    dt = table[,group == colnames(count_table)[i]]
    count_table[,i] = rowSums(dt>0)
  }
  count_table
}

group_sum<- function(table,group,margin=1){
  ### we can use margin to cluster samples(row) or otu(col) according to a grouping factor
  group = as.factor(group)
  level= unique(group)
  if(margin==1){
    dt = sapply(level,function(one_level) rowSums(table[,group==one_level]))
    colnames(dt)=unique(group)
  }else {
    dt = sapply(level,function(one_level) colSums(table[group==one_level,])) %>% t(.)
    rownames(dt)=unique(group)
  }
  as.data.frame(dt)
}

group_mean<- function(table,group,by_row=FALSE){
  ### we can use margin to cluster samples(row) or otu(col) according to a grouping factor
  group = as.factor(group)
  level= unique(group)
  if(by_row==F){
    dt = sapply(level,function(one_level) rowMeans(table[,group==one_level]))
    colnames(dt)=unique(group)
  }else {
    dt = sapply(level,function(one_level) colMeans(table[group==one_level,])) %>% t(.)
    rownames(dt)=unique(group)
  }
  as.data.frame(dt)
}

group_filter = function(table,group,freq=1,cut=1,group_cut=F){
  ### this function is used to filter row data when you group count data table to meet requirement of group count frequency
  ### such as: group_filter(OTU,treatment,freq=3,cut=1,group_cut=F) the otu occured at least 3 times in a level of treatment will be filtered into new datatable.
  ### if set group_cut=T, OTU which occurred lower than cut value in all groups will be discarded
  group= group %>% as.factor(.) %>% droplevels(.) # make factor
  table[is.na(table)]=0 ### must rm NA. since NA will make following judgement NA
  count_table = matrix(ncol= nlevels(group),nrow = nrow(table),dimnames = list(rownames(table),levels(group))) # make count table for grouping
  for (i in 1:ncol(count_table)) {
    dt = table[,group == colnames(count_table)[i]] # make sub-dt for each group
    count_table[,i] = rowSums(dt>0) # sum up this sub-dt at each row from this assigned group
  }
  ### get filter index ###
  if (group_cut==FALSE) {
    otu_check_list = rowSum(count_table) >= freq
    table = table[otu_check_list,]
  }else{
    otu_check_list = rowSum(count_table>=freq) >= cut
    table = table[otu_check_list,]
  }
  ###
  table
}

###### get three column from dist martix
col3m<-function(m,dist_name = 'dist',diag=F){
  ###
  m = as.matrix(m)
  pair = data.frame(row=rownames(m)[row(m)[upper.tri(m,diag = diag)]], 
                    col=colnames(m)[col(m)[upper.tri(m,diag = diag)]], 
                    dist = m[upper.tri(m)])
  colnames(pair)[3]=dist_name
  pair
}
######
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

###### change three column back diag matrix or symmetric matrix 
reshape_col3m<-function(m,name1='name1',name2='name2',value,diag=T,symmetric=F){
  ### this function is very useful to reshpae the col3m data back to a diag matrix or symmetric matrix, which is able to generate dist or run mantel in a proper format
  # m must be dataframe containing pairwise information
  # name1: str, pairwise name1
  # name1: str, pairwise name2
  # value: str, the target index to reshape back
  # format the col3m data
  m = m[,c(name1,name2,value)]
  colnames(m) = c('name1','name2','value')
  # reshape
  reshape.matrix = tidyr::spread(m,key='name2',value='value') %>% as.data.frame(.)
  m.dist =  reshape.matrix[,-1]
  rownames(m.dist) = reshape.matrix[['name1']]
  
  if(diag){
  # add diagonal
  m.dist[,rownames(m.dist)[!rownames(m.dist)%in%colnames(m.dist)]]=NA
  m.dist[colnames(m.dist)[!colnames(m.dist)%in%rownames(m.dist)],]=NA
  m.dist=m.dist[order(rownames(m.dist)),order(colnames(m.dist))]  
  }

  if(symmetric){
  # make symmetric
  m.dist[upper.tri(m.dist)] <- t(m.dist)[upper.tri(m.dist)] 
  }
  as.matrix(m.dist)
}

###### calculate the geospatial distance by lat and lon; need package: dplyr, geoshpere
geo_dist = function(dis,Lon='Lon',Lat='Lat',unit = 'km'){
  dis = dis[,c(Lon,Lat)] %>% as.matrix(.)
  dis = apply(dis, 1, function(y) apply(dis, 1, function(x) distHaversine(y,x)))
  if (unit == 'km'){dis/1000}
}

############ stat tools #############################
######### 
get_aov_cof = function(formula,dt=NULL,cof='Pr(>F)',var=1){
  # var: string or number, to locate the cof of var in the dataframe
  # y : numeric vector
  # x : treat
  # waiting for extend to be complex formula
  result = aov(as.formula(formula),data = dt) %>% summary(.)
  result[[1]][var,cof]
}

anova.adj.R2 <- function(result){
  MSE = tail(result$`Sum Sq`,n=1)/(sum(result$Df)-sum(result$Df[-length(result$Df)])) 
  MST = sum(result$`Sum Sq`)/sum(result$Df)
  1-MSE/MST
}

lmp <- function (modelobject) {
  #if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

lm_par = function(Y,X){
  decay <- lm(Y~X)
  lm.slope = decay$coefficients ['X'] #%>% round(.,3)
  lm.intercept = decay$coefficients ['(Intercept)'] #%>% round(.,3)
  lm.P = lmp(decay) #%>% labelpstar()
  R2= summary(decay)$adj.r.squared %>% round(.,3)
  N = length(summary(decay)$residuals)
  return(c(lm.slope,lm.P,lm.intercept,R2,N))
}

normalization<-function(y){
  x<-y[!is.na(y)]
  x<-(x - min(x)) / (max(x) - min(x))
  y[!is.na(y)]<-x
  return(y)}

#################
get_var_name <- function(var) {
  deparse(substitute(var))
}

is.outlier <- function(x, thres = 3, na.rm = TRUE) {
  abs(x - mean(x, na.rm = na.rm)) > thres * sd(x, na.rm = na.rm)
}

labelp <- function(x){
  if (x <= 0.05){x = "p < 0.05"} 
  else if (x <= 0.1){x = "0.05 < p < 0.1"} 
  else {x = "p > 0.1"}
  P = factor(x,levels=c('p > 0.1','0.05 < p < 0.1','p < 0.05'))
  P
  }
labelpstar <- function(x){
  if (x <= 0.001){P = "***"}
  else if (x <= 0.01){P = "**"}
  else if (x <= 0.05){P = "*"} 
  else if (x <= 0.1){P = "."}
  else {P = ""}
  P
}

###### label the pair by factor
label_pair = function(dt,meta,row='row',col='col',factor,sample){
  # for diagonal-pair data to add pair tag
  factor_row = sapply(dt[[row]], function(x) meta[[factor]][which(sample==x)])
  factor_col = sapply(dt[[col]], function(x) meta[[factor]][which(sample==x)])
  dt$pair=paste(factor_row,'-',factor_col)
}

############ packup functions #############################
### pack up fucntion to add index ###
index_packup <- function(Project,ENV){
  # alpha diversity index
  ENV$Shannon=diversity(t(otu_table(Project)),index='shannon')
  ENV$Richness=specnumber(t(otu_table(Project)))
  #ENV$Chao1 = estimateR(t(otu_table(Project)))['S.chao1',]
  return(ENV)
}

Min3SampleDepth = function(Project) {
  Dt = otu_table(Project)
  Min3SampleDepth = sort(colSums(Dt), decreasing = F)[1:3]
  Min3SampleDepth
}

phylo_pack <- function(dt,taxa,meta,tree,RmSingleSite=F,ReSampleDepth=F){
  ################# remove singleton #############
  dt = dt[rowSums(dt)>=2,]

  ############### phylo object #############
  OTU <- otu_table(dt,taxa_are_rows = TRUE) %>% .[,order(colnames(.))]
  TAX <- taxa %>% as.matrix(.) %>% tax_table(.) 
  ENV <- meta  %>% sample_data(.) %>% .[order(rownames(.)),]
  Project = phyloseq(OTU,TAX,ENV,tree)
  # sum(rownames(OTU) %in% rownames(TAX))
  # sum(colnames(OTU) %in% rownames(ENV))
  # sum(tree$tip.label %in% rownames(TAX))
  ############### Remove Single Site #############
  if (RmSingleSite) {
    index = rowSums(otu_table(Project)>0)>1
    cat("##### Remove Single Site ########")
    cat("### there are ",sum(!index)," once-detected species ###")
    otu_table(Project) = otu_table(Project)[index,]
  }
  ############### other indexes #############
  ENV = sample_data(Project)
  ENV = index_packup(Project,ENV)
  sample_data(Project)=ENV
  ############ Resample Sequencing Depth #############
  if (ReSampleDepth) {
  MinDepth = Min3SampleDepth(Project)
  cat("##### Resample Sequencing Depth #####\n")
  cat("### The lowest 3 depths of sample are ",MinDepth, " ###\n")
  Project = rarefy_even_depth(Project,sample.size=)
  }
  ###
  return(Project)
}

########### utility ##############
### MNTD function define
MT.mntd = function(N,samp,dis,abundance.weighted){
  # require(bigmemory)
  # options(bigmemory.allow.dimnames=TRUE)
  
  sppInSample <- colnames(samp)[samp[N,]>0]
  if (length(sppInSample) > 1) {
    sample.dis <- dis[sppInSample, sppInSample]
    diag(sample.dis) <- NA
    
    #cat("### Now calculating a dissimilarity sub-matrix with ",format(object.size(sample.dis), units = "auto"),"\n")
    
    if (abundance.weighted) {
      mntds <- apply(sample.dis, 1, min, na.rm = TRUE)
      sample.weights <- samp[i, sppInSample]
      mntd <- weighted.mean(mntds, sample.weights)
    } else {
      mntd <- mean(apply(sample.dis, 1, min, na.rm = TRUE))
    }
  }
  else {
    mntd <- NA
  }
  return(mntd)
}

MT.mpd = function(N,samp,dis,abundance.weighted){
  # require(bigmemory)
  # options(bigmemory.allow.dimnames=TRUE)
  
  sppInSample <- colnames(samp)[samp[N,]>0]
  if (length(sppInSample) > 1) {
    sample.dis <- dis[sppInSample, sppInSample]
    diag(sample.dis) <- NA
    
    #cat("### Now calculating a dissimilarity sub-matrix with ",format(object.size(sample.dis), units = "auto"),"\n")
    
    if (abundance.weighted) {
      sample.weights <- t(as.matrix(samp[i, sppInSample, drop = FALSE])) %*% 
        as.matrix(samp[i, sppInSample, drop = FALSE])
      mpd[i] <- weighted.mean(sample.dis, sample.weights)
    } else {
      mpd <- mean(apply(sample.dis, 1, mean, na.rm = TRUE))
    }
  }
  else {
    mpd <- NA
  }
  return(mpd)
}

MT.mpd.mntd = function(N,samp,pd.wd,pd.desc,pd.spname,abundance.weighted){
  require(bigmemory)
  # re-attch dissimilarity matrix
  dis = attach.big.matrix(dget(paste0(pd.wd, "/",pd.desc)))
  options(bigmemory.allow.dimnames=TRUE)
  dimnames(dis) = list(pd.spname,pd.spname)
  #cat("### Attached a dissimilarity matrix with ",format(object.size(dis), units = "auto"),"\n")
  
  sppInSample <- colnames(samp)[samp[N,]>0]
  if (length(sppInSample) > 1) {
    #cat("# subset matrix start:",format(Sys.time(), "%H:%M:%S"),"\n")
    sample.dis <- dis[sppInSample, sppInSample]
    diag(sample.dis) <- NA
    
    #cat("### Size of sub-matrix: ",format(object.size(sample.dis), units = "auto"),"\n")
    mpd <- apply(sample.dis, 1, mean, na.rm = TRUE) %>% mean(.,na.rm = T)
    mntd <- apply(sample.dis, 1, min, na.rm = TRUE) %>% mean(.,na.rm = T)
    # mntd = unique(mntd)
    #cat("# subset matrix end:",format(Sys.time(), "%H:%M:%S"),"\n")
  }  else {
    mntd = mpd = NA
  }
  rm(sample.dis)
  gc()
  return(c(mpd,mntd))
}

