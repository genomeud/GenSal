# Run with --help flag for help.
# Modified 08/01/2020 by Fabio Marroni
#We compared atlantica and marmorata transcripts. So the transcript "atlantica private" are those absent in marmorata.

suppressPackageStartupMessages({
  library(optparse)
})


option_list = list(
  make_option(c("-I", "--infile"), type="character", default="",
              help="HapFLK output file", metavar="character"),
  make_option(c("-P", "--pvalue"), type="numeric", default=0.05, 
              help="pvalue to be considered significant [default= %default]", metavar="character"),
  make_option(c("-D", "--maxdist"), type="numeric", default=1000000, 
              help="Maximum distance for two signal to be clustered [default= %default]", metavar="character"),
  make_option(c("-S", "--minsnp"), type="numeric", default=2, 
              help="Minimum number of consecutive significant SNPs to be included in results [default= %default]", metavar="character"),
  make_option(c("-L", "--maxclusters"), type="numeric", default=50, 
              help="Maximum number of clusters to be recorded for comparison [default= %default]", metavar="character"),
  make_option(c("-N", "--numsim"), type="numeric", default=100, 
              help="Number of simulations to run [default= %default]", metavar="character"),
  make_option(c("-O", "--outfile"), type="character", default="", 
              help="output file (containing top results of HapFLK and formatted for 07_KEGG_enrichment.r) [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

  if (is.null(opt$infile)) {
  stop("WARNING: No infile specified with '-I' flag.")
} else {  cat ("infile is ", opt$infile, "\n")
  infile <- opt$infile  
  }

  if (is.null(opt$pvalue)) {
  stop("WARNING: No pvalue specified with '-P' flag.")
} else {  cat ("pvalue is ", opt$pvalue, "\n")
  pvalue <- opt$pvalue  
  }

  if (is.null(opt$maxdist)) {
  stop("WARNING: No maxdist specified with '-D' flag.")
} else {  cat ("maxdist is ", opt$maxdist, "\n")
  maxdist <- opt$maxdist  
  }

  if (is.null(opt$minsnp)) {
  stop("WARNING: No minsnp specified with '-S' flag.")
} else {  cat ("minsnp is ", opt$minsnp, "\n")
  minsnp <- opt$minsnp  
  }

  if (is.null(opt$maxclusters)) {
  stop("WARNING: No maxclusters specified with '-L' flag.")
} else {  cat ("maxclusters is ", opt$maxclusters, "\n")
  maxclusters <- opt$maxclusters  
  }

  if (is.null(opt$numsim)) {
  stop("WARNING: No numsim specified with '-N' flag.")
} else {  cat ("numsim is ", opt$numsim, "\n")
  numsim <- opt$numsim  
  }

  if (is.null(opt$outfile)) {
  stop("WARNING: No outfile specified with '-O' flag.")
} else {  cat ("outfile is ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }

extract_HFLK_sig<-function(infile,threshold,outfile,pvalue,maxdist,minsnp,maxclusters,numsim)
{
library("data.table")
hFLK<-fread(infile,data.table=F)
#Keep a copy of all hFLK values, we will need to compute some statistics
allhFLK<-hFLK
#Remove one useles column
hFLK$rs<-hFLK$hapflk<-hFLK$hapflk_scaled<-NULL
hFLK$chr<-gsub("b'","",hFLK$chr)
hFLK$chr<-gsub("'","",hFLK$chr)
#Just to be sure, sort by chromosome and position
hFLK<-hFLK[order(hFLK$chr,hFLK$pos),]
#Group regions in a similar way to what we did for ZHp.
#Here we do not have a "real" explicit window.
#We therefore take note of consecutive significant pvalues and call a window as "signficant" if it contains more than X (to decide)
#Significant pvalues not interrupted by any not signficant pvalue.
#This still has to be implemented.
#I think it is ok to collapse different GBS loci if they show a consistent signal, because they might belong to the same haplotype.

hFLK$group<-NA
mygroup<-1
#Ouch another horrible loop!
#Here, I just assign positions to groups.
#Every time the the chromosome changes, the distance with the previous result is greater than 1Mb (default, will be a parameter),
# or the result is not significant, the group changes.

clusters<-matrix(0,nrow=numsim,ncol=maxclusters)
oriFLK<-hFLK
for(bbb in 1:nrow(clusters))
{
cat("Simulation",bbb,"of",nrow(clusters),"\n")
hFLK<-oriFLK
if(bbb!=1) hFLK$pvalue<-sample(hFLK$pvalue)
for(aaa in 2:nrow(hFLK))
{
	if(hFLK$chr[aaa-1]!=hFLK$chr[aaa]) 
	{
		mygroup<-mygroup+1
		if(hFLK$pvalue[aaa]<=pvalue) hFLK$group[aaa]<-mygroup
		next
	}
	if(hFLK$pos[aaa]-hFLK$pos[aaa-1]>maxdist)
	{
		mygroup<-mygroup+1
		if(hFLK$pvalue[aaa]<=pvalue) hFLK$group[aaa]<-mygroup
		next
	}
	if(hFLK$pvalue[aaa]>pvalue)
	{
		mygroup<-mygroup+1
		next
	}
	hFLK$group[aaa]<-mygroup
}

#We now remove all non-significant results, i.e. those with NA in the group column
hFLK<-hFLK[!is.na(hFLK$group),]
#Calculate start, end, and length of windows in bp
ll<-aggregate(hFLK$pos,by=list(hFLK$group),FUN="range")
ll$length<-abs(ll$x[,1]-ll$x[,2])
ll$start<-ll$x[,1]
ll$end<-ll$x[,2]
ll<-ll[,c("Group.1","start","end","length")]
finalFLK<-aggregate(hFLK[,!names(hFLK)%in%"chr"],by=list(hFLK$group),FUN="median")
finalFLK<-merge(ll,finalFLK)
#Calculate length of windows in number of SNPs
ll<-aggregate(hFLK$pos,by=list(hFLK$group),FUN="length")
finalFLK<-merge(ll,finalFLK)
setnames(finalFLK,"x","nSNP")
#Record the top signal (best pvalue) for each window
ll<-aggregate(hFLK$pvalue,by=list(hFLK$group),FUN="min")
finalFLK<-merge(ll,finalFLK)
setnames(finalFLK,"x","top")
#Reassign chromosomes, based on group
temp<-hFLK[,c("chr","group")]
temp<-temp[!duplicated(temp),]
finalFLK<-merge(temp,finalFLK,by="group")
finalFLK$log10p<-(-log10(finalFLK$pvalue))
finalFLK$top<-(-log10(finalFLK$top))
finalFLK$group<-finalFLK$pvalue<-finalFLK$Group.1<-NULL
#Select only useful columns
finalFLK<-finalFLK[,c("chr","start","end","length","nSNP","log10p","top")]
#Order
finalFLK<-finalFLK[order(finalFLK$top,decreasing=T),]
#Remove windows composed by less than minsnp SNPs
finalFLK<-finalFLK[finalFLK$nSNP>=minsnp,]
endfill<-min(nrow(finalFLK),maxclusters)
clusters[bbb,1:endfill]<-sort(finalFLK$nSNP,decreasing=T)[1:endfill]
}

sumcl<-clusters[1:2,]
for(ddd in 1:ncol(sumcl))
{
sumcl[2,ddd]<-sum(clusters[2:nrow(clusters),]>=sumcl[1,ddd])/(nrow(clusters)*ncol(clusters))
}
fincl<-data.frame(Size=sumcl[1,],pvalue=sumcl[2,])
browser()
write.table(fincl,outfile,quote=F,sep="\t",row.names=F)

}

extract_HFLK_sig(infile=infile,pvalue=pvalue,maxdist=maxdist,minsnp=minsnp,outfile=outfile,maxclusters=maxclusters,numsim=numsim)
