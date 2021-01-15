# Copyright:	    Gabriele Magris & Fabio Marroni 2020
# Aim:              Draw treemix plot 
# To add:		     
# Suggestions: 
# Fixes:  

suppressPackageStartupMessages(library("optparse"))

option_list = list(
  make_option(c("-d", "--wdir"), action="store", default=NULL, type='character', 
              dest="wd", help="Directory where file are located  [%default]"),
  make_option(c("-f", "--funcdir"), action="store", default=NULL, type='character', 
              dest="funcdir", help="Directory where functions are located  [%default]"),
  make_option(c("-p", "--prefix"), action="store", default=NULL, type='character',
              help="Prefix of input file  [%default]")
              )
opt = parse_args(OptionParser(option_list=option_list,description="\nCompute SNP phasing from re-sequenced selfed individuals - only heterozygous regions"))
# if(length(opt$args) != 1) {
	# print_help(opt)
# }
# print(str(opt))

wd=opt$wd
funcdir=opt$funcdir
prefix=opt$prefix


# draw plot 
draw_treemix_plot<-function(wd="",prefix="",funcdir="")
{
    setwd(wd)
    source(paste(funcdir,"/a04_treemix_sub.r",sep=""))
    dir.create("../plots/",recursive=T,showWarnings=F)
    infile<-prefix
    
    # set name to save plot 
    out_plot<-paste("../plots/",infile,".pdf",sep="")
    pdf(out_plot,width=6,height=6)
    plot_tree(infile,disp=0.0003,mbar=F,cex=0.6)
    
    plot_resid(infile,"../groups/poporder_K5_ddrad.txt",cex=0.6)
    dev.off()
 
}

draw_treemix_plot(wd=wd,prefix=prefix,funcdir=funcdir)