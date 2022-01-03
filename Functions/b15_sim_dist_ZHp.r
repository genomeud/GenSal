# Run with --help flag for help.
# Modified 08/01/2020 by Fabio Marroni

suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-Z", "--ZHpfile"), type="character", default="",
              help="Output file created by ZHp, to be used as input of this function", metavar="character"),
  make_option(c("-T", "--threshold"), type="numeric", default=2.81, 
              help="Z-score significance threshold [default= %default]", metavar="character"),
  make_option(c("-L", "--maxlength"), type="numeric", default=27, 
              help="Max number of windows to investigate [default= %default]", metavar="character"),
  make_option(c("-S", "--nsim"), type="numeric", default=1000, 
              help="Number of simulations [default= %default]", metavar="character"),
  make_option(c("-O", "--outfile"), type="character", default="", 
              help="output file name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

  if (is.null(opt$ZHpfile)) {
  stop("WARNING: No ZHpfile specified with '-Z' flag.")
} else {  cat ("ZHpfile is ", opt$ZHpfile, "\n")
  ZHpfile <- opt$ZHpfile  
  }

  if (is.null(opt$threshold)) {
  stop("WARNING: No threshold specified with '-T' flag.")
} else {  cat ("threshold is ", opt$threshold, "\n")
  threshold <- opt$threshold  
  }

  if (is.null(opt$maxlength)) {
  stop("WARNING: No maxlength specified with '-L' flag.")
} else {  cat ("maxlength is ", opt$maxlength, "\n")
  maxlength <- opt$maxlength  
  }

  if (is.null(opt$nsim)) {
  stop("WARNING: No nsim specified with '-S' flag.")
} else {  cat ("nsim is ", opt$nsim, "\n")
  nsim <- opt$nsim  
  }

  if (is.null(opt$outfile)) {
  stop("WARNING: No outfile specified with '-O' flag.")
} else {  cat ("outfile is ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }

sim_dist_ZHp<-function(ZHpfile,threshold,nsim,maxlength,outfile)
{
library("data.table")
nsim<-nsim+1
clusters<-matrix(0,nrow=nsim,ncol=maxlength)
ZHp<-fread(ZHpfile,data.table=F)
#Remove some lines with NAs value in ZHp
ZHp<-na.omit(ZHp)
#Keep a copy of all ZHp values, we will need to compute some statistics
allZHp<-ZHp
for(bbb in 1:nsim)
{
ZHp<-allZHp
if(bbb>1) ZHp$ZHp<-sample(ZHp$ZHp)
#Only select values that exceed the threshold
ZHp<-ZHp[abs(ZHp$ZHp)>=threshold,]
#We count the number of overlapping windows in which at results exced the threshold
ZHp$group<-NA
mygroup<-1
for(aaa in 2:nrow(ZHp))
{
if(ZHp$chr[aaa-1]!=ZHp$chr[aaa]) 
{
mygroup<-mygroup+1
next
}
temp.ov<-count.overlap(ZHp$start[aaa-1],ZHp$end[aaa-1],ZHp$start[aaa],ZHp$end[aaa])
#If two consecutive windows overlap and have same ZHp sign we assign them to the same group
if(temp.ov>0 & (ZHp$ZHp[aaa-1]*ZHp$ZHp[aaa])>0) ZHp$group[aaa-1]<-ZHp$group[aaa]<-mygroup
if(temp.ov<1) mygroup<-mygroup+1
}
#Here we summarize (median) ZHp where ZHp exceeds threshold in at least 2 overlapping windows 
ZHp<-na.omit(ZHp)
if (nrow(ZHp)<1) next
aggname<-aggregate(ZHp$chr,by=list(ZHp$group),FUN=unique)
setnames(aggname,"x","chr")
aggstart<-aggregate(ZHp$start,by=list(ZHp$group),FUN=min)
setnames(aggstart,"x","start")
aggend<-aggregate(ZHp$end,by=list(ZHp$group),FUN=max)
setnames(aggend,"x","end")
#This would compute only the median value of ZHp across significant windows
# aggzhp<-aggregate(ZHp$ZHp,by=list(ZHp$group),FUN=median)
# setnames(aggzhp,"x","ZHp")
agglength<-aggregate(ZHp$ZHp,by=list(ZHp$group),FUN=length)
setnames(agglength,"x","length")

mysummary<-merge(aggname,agglength,sort=F)
mysummary<-merge(mysummary,aggstart,sort=F)
mysummary<-merge(mysummary,aggend,sort=F)

#mysummary<-merge(mysummary,aggzhp,sort=F)
#Here we use allZHp to compute the median value of ZHp across all the windows included in the selected windows (i.e. we also include non signficant reults)
mysummary$ZHp<-NA
for(ccc in 1:nrow(mysummary))
{
pp<-allZHp[(allZHp$chr==mysummary$chr[ccc])&(allZHp$start>=mysummary$start[ccc])&(allZHp$end<=mysummary$end[ccc]),]
mysummary$ZHp[ccc]<-median(pp$ZHp)
}

endfill<-min(nrow(ZHp),maxlength)
clusters[bbb,1:endfill]<-sort(mysummary$length,decreasing=T)[1:endfill]
}
clusters[is.na(clusters)]<-0
sumcl<-clusters[1:2,]
for(ddd in 1:ncol(sumcl))
{
sumcl[2,ddd]<-sum(clusters[2:nrow(clusters),]>=sumcl[1,ddd])/(nrow(clusters)*ncol(clusters))
}
fincl<-data.frame(Size=sumcl[1,],pvalue=sumcl[2,])

write.table(fincl,outfile,quote=F,sep="\t",row.names=F)
}




#Here I just write all the instances belonging to 


count.overlap<-function(start.1=1,end.1=4845859,start.2=1,end.2=84945)
{
	max.start<-max(start.1,start.2)
	min.end<-min(end.1,end.2)
	overlap<-max(0,(1+min.end-max.start))
	overlap
}

count.overlap.vector<-function(x)
{
	max.start<-max(x[1],x[3])
	min.end<-min(x[2],x[4])
	overlap<-max(0,(1+min.end-max.start))
	overlap
}


sim_dist_ZHp(ZHpfile=ZHpfile,threshold=threshold,outfile=outfile,maxlength=maxlength,nsim=nsim)
