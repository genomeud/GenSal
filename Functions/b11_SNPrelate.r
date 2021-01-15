# Run with --help or -h flag for help.
# Written 04/07/2020 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-G", "--gdsfile"), type="character", default="",
		help="SNPRelate gds file, to be used as input [default= %default]", metavar="character"), 
  make_option(c("-P", "--popfile"), type="character", default="",
		help="File containing population assignment [default= %default]", metavar="character"), 
  make_option(c("-B", "--between"), type="character", default="",
		help="Between population IBD graph file [default= %default]", metavar="character"), 
  make_option(c("-W", "--within"), type="character", default="",
		help="Within population IBD graph file [default= %default]", metavar="character"), 
  make_option(c("-H", "--highcovdir"), type="character", default="",
		help="Folder containing only samples with cov > 5 (only needed to filter out low cov) [default= %default]", metavar="character"), 
  make_option(c("-M", "--maf"), type="numeric", default=0.05,
		help="Minimum MAF to include a SNP in analysis [default= %default]", metavar="character"), 
  make_option(c("-S", "--missingness"), type="numeric", default=0.05,
		help="Maximum missingness to include a SNP in analysisi [default= %default]", metavar="character"), 
  make_option(c("-L", "--LD"), type="numeric", default=0.2,
		help="Threshold for LD pruning [default= %default]", metavar="character"), 
  make_option(c("-I", "--ibdmatfile"), type="character", default="", 
		help="IBD matrix output file (computed using SNPRelate) [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")


if (is.null(opt$gdsfile)) {
  stop("WARNING: No gds file specified with '-G' flag.")
} else {  cat ("gds file ", opt$gdsfile, "\n")
  gdsfile <- opt$gdsfile  
  }

  if (is.null(opt$popfile)) {
  stop("WARNING: No popfile specified with '-P' flag.")
} else {  cat ("popfile is", opt$popfile, "\n")
  popfile <- opt$popfile  
  }

  if (is.null(opt$between)) {
  stop("WARNING: No between file specified with '-B' flag.")
} else {  cat ("between is", opt$between, "\n")
  between <- opt$between  
  }

  if (is.null(opt$within)) {
  stop("WARNING: No within file specified with '-W' flag.")
} else {  cat ("within is", opt$within, "\n")
  within <- opt$within  
  }

  if (is.null(opt$highcovdir)) {
  stop("WARNING: No highcov dir specified with '-H' flag.")
} else {  cat ("highcovdir is", opt$highcovdir, "\n")
  highcovdir <- opt$highcovdir  
  }

  if (is.null(opt$maf)) {
  stop("WARNING: No maf specified with '-M' flag.")
} else {  cat ("maf is", opt$maf, "\n")
  maf <- opt$maf  
  }

  if (is.null(opt$missingness)) {
  stop("WARNING: No missingness specified with '-S' flag.")
} else {  cat ("missingness is", opt$missingness, "\n")
  missingness <- opt$missingness  
  }

  if (is.null(opt$LD)) {
  stop("WARNING: No LD specified with '-L' flag.")
} else {  cat ("LD is", opt$LD, "\n")
  LD <- opt$LD  
  }

  if (is.null(opt$ibdmatfile)) {
  stop("WARNING: No ibdmatfile specified with '-I' flag.")
} else {  cat ("ibdmatfile is", opt$ibdmatfile, "\n")
  ibdmatfile <- opt$ibdmatfile  
  }

run_ibd<-function(gdsfile,coefficient,popfile,between,within,ibdmatfile,maf,missingness,LD,highcovdir)
{
library(data.table)
library(RColorBrewer)
library(gdsfmt)
library(SNPRelate)
#Load gds file (has to be previously saved)
genofile <- snpgdsOpen(gdsfile)
#Keep only high coverage samples. We read the names from an intermediate results directory
keep.me<-gsub(".reference.cov5.info50.no_low_cov.txt","",dir(highcovdir,pattern=".reference.cov5.info50.no_low_cov.txt"))

#Prune for LD (otherwise IBD estimates might be inflated and running time longer)
snpset <- snpgdsLDpruning(genofile, ld.threshold=LD,autosome.only=F)
snpset.id<-unlist(snpset)

sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
#Compute IBD and write IBD matrix just in case...
ibd <- snpgdsIBDMLE(genofile, maf=maf, missing.rate=missingness, sample.id=keep.me,num.thread=8,snp.id=snpset.id,autosome.only=F)
ibd.coeff <- snpgdsIBDSelection(ibd)
ibd.coeff<-ibd.coeff[order(ibd.coeff$kinship,decreasing=T),]

write.table(ibd.coeff,ibdmatfile,row.names=F,quote=F,sep="\t")

#Assign populations to individuals
mypop<-fread(popfile,data.table=F)
mypop$final_population[mypop$final_population=="fario_atlantica"]<-"AT"
mypop$final_population[mypop$final_population=="fario_med_island"]<-"MI"
mypop$final_population[mypop$final_population=="fario_med_peninsula"]<-"MM"
mypop$final_population[mypop$final_population=="carpione"]<-"CA"
mypop$final_population[mypop$final_population=="marmorata"]<-"MA"

ibd.coeff$POP2<-ibd.coeff$POP1<-NULL
ibd.coeff$POP1<-mypop$final_population[match(ibd.coeff$ID1,mypop$sample_ID)]
ibd.coeff$POP2<-mypop$final_population[match(ibd.coeff$ID2,mypop$sample_ID)]

#Only select within population comparisons
intra<-ibd.coeff[ibd.coeff$POP1==ibd.coeff$POP2,]
#Plot within pop relatedness
intra$POP1<-factor(intra$POP1,levels=c("MI","AT","MM","MA","CA"))
#intra$ind1.id<-factor(intra$ind1.id,levels=c("Mediterranea\nIsland","Atlantica","Mediterranea\nMainland","Marmoratus","Carpione"))
color<-c("blue","green","orangered","gray68","orchid2")
pdf(within,height=6)
par(mgp=c(3,1.2,0))
boxplot(intra$kinship~intra$POP1,col=color,cex.axis=0.8,ylab="IBD",xlab="Sample")
dev.off()

#Only select between population comparisons
inter<-ibd.coeff[ibd.coeff$POP1!=ibd.coeff$POP2,]
inter$POP1<-factor(inter$POP1,levels=c("MI","AT","MM","MA","CA"))
inter<-inter[order(inter$POP1),]
inter$pair<-paste(as.character(inter$POP1),inter$POP2,sep="\n")
color<-colorRampPalette(brewer.pal(9,"Blues"))(length(unique(inter$pair)))
pdf(between,height=6)
par(mgp=c(3,1.2,0))
boxplot(inter$kinship~inter$pair,col=color,cex.axis=0.8,ylab="IBD",xlab="Comparison")
dev.off()

}
run_ibd(gdsfile=gdsfile,coefficient=coefficient,popfile=popfile,between=between,within=within,ibdmatfile=ibdmatfile,maf=maf,missingness=missingness,LD=LD,highcovdir=highcovdir)
