source('0_function.R')
library(iCAMP)
library(NST)
library(ieggr)
library(tidyr)
cat("##############################\n")


load('7_icamp_ITS_ASV/icamp_preload_ITS_ASV.rda')
#################### save dir #######################
folder = paste0('7_icamp_',prefix)
dir.create(folder)
###############################################
##############################
####
# computational setting
nworker   =28
memory.G  =120
rand.time = 100
ds = 0.2 # significant distance threshold
bin.size.limit = 48
abcut=5 # you may remove some species, if they are too rare to perform reliable correlation test. ###
omit.option = "omit"
sig.index="Confidence"

cat("### threads:", nworker,"\n")
cat("### memory:", memory.G,"\n")
cat("### bin.size.limit:", bin.size.limit,"\n")
cat("### omit.option:", omit.option,"\n")

###############################################
#################### save dir #######################
folder = paste0('7_icamp_',prefix)
dir.create(folder)

intermediate_file = paste0(folder,'/intermediate_file')
unlink(intermediate_file, recursive = T, force = FALSE)
dir.create(intermediate_file)

note ='test' # for optional subset
##################################################
cat("### phylo distance matrix constructing\n")

pd.big=pdist.big(tree = Tre, wd=intermediate_file, nworker = nworker, memory.G = memory.G)
save(list = c('pd.big'),file = paste0(intermediate_file,'/pd.big'))

# phylobin
cat("### phylo binning\n")
phylobin=taxa.binphy.big(tree = Tre, pd.desc = pd.big$pd.file,pd.spname = pd.big$tip.label,
                         pd.wd = pd.big$pd.wd, ds = ds, bin.size.limit = bin.size.limit,
                         nworker = nworker)
save(file = paste0(intermediate_file,'/phylobin'),list = c('phylobin'))

# bin parameters or vectors
cat("### phylo bin signals calculating by mantel\n")
niche.dif=dniche(env = Var,comm = t(OTU),method = "niche.value",
                 nworker = nworker,out.dist=F,bigmemo=TRUE,
                 nd.wd=intermediate_file)

sp.bin=phylobin$sp.bin[,3,drop=FALSE]
sp.ra=colMeans(t(OTU)/rowSums(t(OTU)))
commc=t(OTU)[,colSums(t(OTU))>=abcut,drop=FALSE]
dim(commc)
spname.use=colnames(commc)
binps=ps.bin(sp.bin = sp.bin,
             sp.ra = sp.ra,
             spname.use = spname.use,
             pd.desc = pd.big$pd.file, 
             pd.spname = pd.big$tip.label, 
             pd.wd = pd.big$pd.wd,
             nd.list = niche.dif$nd,
             nd.spname = niche.dif$names,
             ndbig.wd = niche.dif$nd.wd,
             cor.method = "spearman",
             r.cut = 0.2, 
             p.cut = 0.05, 
             min.spn = bin.size.limit)

### write bin signal result
if(file.exists(paste0(folder,'/',prefix,".PhyloSignalSummary.csv"))){appendy=TRUE;col.namesy=FALSE}else{appendy=FALSE;col.namesy=TRUE}

write.table(data.frame(ds=ds,n.min=bin.size.limit,binps$Index),file = paste0(folder,'/',prefix,".PhyloSignalSummary.csv"),
            append = appendy, quote=FALSE, sep=",", row.names = FALSE,col.names = col.namesy)

if(file.exists(paste0(folder,'/',prefix,".PhyloSignalDetail.csv"))){appendy2=TRUE;col.namesy2=FALSE}else{appendy2=FALSE;col.namesy2=TRUE}

write.table(data.frame(ds=ds,n.min=bin.size.limit,binID=rownames(binps$detail),binps$detail),file = paste0(folder,'/',prefix,".PhyloSignalDetail.csv"),
            append = appendy2, quote = FALSE, sep = ",", row.names = FALSE, col.names = col.namesy2)

########## bin camp
cat("### phylo bin icamps calculating\n")
comm=t(OTU)
icres=icamp.big(comm=comm, 
                pd.desc = pd.big$pd.file, 
                pd.spname=pd.big$tip.label,
                pd.wd = pd.big$pd.wd, 
                rand = rand.time, 
                tree=Tre,
                prefix = prefix, 
                ds = 0.2, 
                pd.cut = NA, 
                sp.check = TRUE,
                phylo.rand.scale = "within.bin", 
                taxa.rand.scale = "across.all",
                phylo.metric = "bMPD", 
                sig.index=sig.index, 
                bin.size.limit = bin.size.limit, 
                nworker = nworker, 
                memory.G = memory.G, 
                rtree.save = FALSE, 
                detail.save = TRUE, 
                qp.save = T, 
                detail.null = T, 
                ignore.zero = TRUE, 
                output.wd = intermediate_file, 
                correct.special = TRUE, 
                unit.sum = rowSums(t(OTU)), 
                special.method = "depend",
                ses.cut = 1.96, 
                rc.cut = 0.95, 
                conf.cut=0.975, 
                omit.option = "no",
                meta.ab = NULL)

save(file = paste0(intermediate_file,'/','icres'),list = c('icres'))


# iCAMP bin level statistics
cat("### iCAMP bin level statistics\n")
icbin=icamp.bins(icamp.detail = icres$detail,
                 treat = treat,
                 clas=TAX,
                 silent=FALSE, 
                 boot = TRUE,
                 rand.time = rand.time,
                 between.group = TRUE)
# Bootstrapping test
cat("### phylo bin icamps bootstrapping\n")

icamp.result=icres$CbMPDiCBraya
icboot=iCAMP::icamp.boot(icamp.result = icamp.result,
                         treat = treat,
                         rand.time = rand.time,
                         compare = TRUE,
                         silent = FALSE,
                         between.group = TRUE,
                         ST.estimation = TRUE)

# output files:
# Test.iCAMP.Summary.rda: the object "icbin" saved in R data format. see help document of the function icamp.bins for description of each element in the object.
# Test.ProcessImportance_EachGroup.csv: Relative importance of each process in governing the turnovers in a group of samples.
# Test.ProcessImportance_EachBin_EachGroup.csv: Relative importance of each process in governing the turnovers of each bin among a group of samples.
# Test.ProcessImportance_EachTurnover.csv: Relative importance of each process in governing the turnovers between each pair of communities (samples).
# Test.BinContributeToProcess_EachGroup.csv: Bin contribution to each process, measuring the contribution of each bin to the relative importance of each process in the assembly of a group of communities.
# Test.Taxon_Bin.csv: a matrix showing the bin ID and classification information for each taxon.
# Test.Bin_TopTaxon.csv: a matrix showing the bin relative abundance; the top taxon ID, percentage in bin, and classification; the most abundant name at each phylogeny level in the bin.
save(icbin,file = paste0(folder,'/',prefix,".iCAMP.Summary.rda")) # just to archive the result. rda file is automatically compressed, and easy to load into R.
write.csv(icbin$Pt,file = paste0(folder,'/',prefix,".ProcessImportance_EachGroup.csv"),row.names = FALSE)
write.csv(icbin$Ptk,file = paste0(folder,'/',prefix,".ProcessImportance_EachBin_EachGroup.csv"),row.names = FALSE)
write.csv(icbin$Ptuv,file = paste0(folder,'/',prefix,".ProcessImportance_EachTurnover.csv"),row.names = FALSE)
write.csv(icbin$BPtk,file = paste0(folder,'/',prefix,".BinContributeToProcess_EachGroup.csv"),row.names = FALSE)
write.csv(data.frame(ID=rownames(icbin$Class.Bin),icbin$Class.Bin,stringsAsFactors = FALSE),
          file = paste0(folder,'/',prefix,".Taxon_Bin.csv"),row.names = FALSE)
write.csv(icbin$Bin.TopClass,file = paste0(folder,'/',prefix,".Bin_TopTaxon.csv"),row.names = FALSE)

save(icboot,file=paste0(folder,'/',prefix,".iCAMP.Boot.",note,".rda"))
write.csv(icboot$summary,file = paste0(folder,'/',prefix,".iCAMP.BootSummary.",note,".csv"),row.names = FALSE)
write.csv(icboot$compare,file = paste0(folder,'/',prefix,".iCAMP.Compare.",note,".csv"),row.names = FALSE)

rm(list=ls())

