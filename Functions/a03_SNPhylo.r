# Copyright:	    Gabriele Magris & Fabio Marroni 2020
# Aim:              Compute bootstrap analysis on phylo tree
# To add:		     
# Suggestions: 
# Fixes:  

suppressPackageStartupMessages(library("optparse"))

option_list = list(
  make_option(c("-i", "--inpath"), action="store", default=NULL, type='character', 
              dest="wd", help="Directory where file are located  [%default]"),
  make_option(c("-p", "--prefix"), action="store", default="salmo_trutta_id1492.cov5.info50.no_low_cov", type='character',
              dest="prefix",help="File name prefix [%default]")
			 )

opt_parser = OptionParser(option_list=option_list,description="\nCompute bootstrap analysis on phylogenetic tree")
opt = parse_args(OptionParser(option_list=option_list,description="\nCompute bootstrap analysis on phylogenetic tree"))

if(is.null(opt$wd)) {
	print_help(opt_parser)
    stop("At least one argument must be supplied", call.=FALSE)
}

inpath=opt$wd
file.prefix=opt$prefix

perform_bootstrap_on_phylo_tree<-function(inpath,file.prefix)
{
# load the necessary library
library(phangorn)
library(methods)
bs_image_file.name <- paste(file.prefix, ".bs.png", sep = "")
bs_tree_file.name  <- paste(file.prefix, ".bs.tree", sep = "")

phylip <- read.phyDat(paste(inpath,"/",file.prefix,".phylip.txt",sep=""), format="phylip", type="DNA")
newick <- read.tree(paste(inpath,"/",file.prefix,".ml.bk.tree",sep=""))
fit <- pml(newick, phylip)
num.bs_sample <- 100

set.seed(1)
bs <- bootstrap.pml(fit, bs = num.bs_sample, optNni=TRUE, multicore=TRUE,mc.cores=1)
png(filename = paste(inpath,"/",bs_image_file.name,sep=""), width = 1000, height = 1000,type="cairo")
options(warn=-1)
bs_tree <- plotBS(fit$tree, bs, cex = 1, edge.width = 2)
options(warn=0)
dev.off()
write.tree(bs_tree, file = paste(inpath,"/",bs_tree_file.name,sep=""))
}


perform_bootstrap_on_phylo_tree(inpath=inpath,file.prefix=file.prefix)

