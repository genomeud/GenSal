# Run with --help flag for help.
# Modified 08/01/2020 by Fabio Marroni

suppressPackageStartupMessages({
  library(optparse)
})


option_list = list(
  make_option(c("-I", "--indir"), type="character", default="",
              help="Directory containing ABBABABA output files to be sued as input", metavar="character"),
  make_option(c("-E", "--estimator"), type="character", default="D", 
              help="Estimator for which we need the summary. Possible choices, D and f [default= %default]", metavar="character"),
  make_option(c("-R", "--rule"), type="character", default="non-overlapping", 
              help="Rule to select windows on which to run the jacknife. Currently only the default available [default= %default]", metavar="character"),
  make_option(c("-O", "--outfile"), type="character", default="", 
              help="output file name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

  if (is.null(opt$indir)) {
  stop("WARNING: No indir specified with '-I' flag.")
} else {  cat ("indir is ", opt$indir, "\n")
  indir <- opt$indir  
  }

  if (is.null(opt$rule)) {
  stop("WARNING: No rule specified with '-R' flag.")
} else {  cat ("rule is ", opt$rule, "\n")
  rule <- opt$rule  
  }

  if (is.null(opt$estimator)) {
  stop("WARNING: No estimator specified with '-E' flag.")
} else {  cat ("estimator is ", opt$estimator, "\n")
  estimator <- opt$estimator  
  }

  if (is.null(opt$outfile)) {
  stop("WARNING: No outfile specified with '-O' flag.")
} else {  cat ("outfile is ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }

run_jack<-function(indir,rule,estimator,outfile)
{
library("data.table")
library(bootstrap)
fullfile<-dir(indir,full.names=T)
mysummary<-data.frame(Sample=basename(fullfile),Est=NA,LCL=NA,UCL=NA,jack=NA,jack.se=NA,stringsAsFactors=F)
for(ccc in 1:length(fullfile))
{
	ABBA<-fread(fullfile[ccc],data.table=F)
	#Select the chosen estimator
	ifelse(estimator=="D",ABBA$Est<-ABBA$D,ABBA$Est<-ABBA$fd)
	#We remove non-informative windows, we don't want to select those!
	ABBA<-ABBA[!is.na(ABBA$Est),]
	ABBA$keep<-0
	ABBA$keep[1]<-1
	last<-1
	if(rule!="non-overlapping") stop("Can you read the help text? Only 'non-overlapping' method is currently available")
	for (aaa in 2:nrow(ABBA))
	{
		#If the window we are testing is on a different chromosome or doesn't overlap the last selected window
		#We mark it as to keep and use it for the following tests
		if((ABBA$scaffold[last]!=ABBA$scaffold[aaa])|(count.overlap(ABBA$start[last],ABBA$end[last],ABBA$start[aaa],ABBA$end[aaa])==0))
		{
			ABBA$keep[aaa]<-1+ABBA$keep[last]
			last<-aaa
		}
	}
	ABBA<-ABBA[ABBA$keep>0,]
	ABBA$keep<-NULL

	#I define my confidence intervals for the estimator
rr<-bootstrap(ABBA$Est,1000,theta=median)
pp<-jackknife(ABBA$Est,theta=median)
mysummary$Est[ccc]<-median(rr$thetastar)
mysummary$LCL[ccc]<-quantile(rr$thetastar,0.025)
mysummary$UCL[ccc]<-quantile(rr$thetastar,0.975)
mysummary$jack[ccc]<-median(pp$jack.values)
mysummary$jack.se[ccc]<-pp$jack.se
}
write.table(mysummary,outfile,quote=F,sep="\t",row.names=F)
}




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


run_jack(indir=indir,rule=rule,estimator=estimator,outfile=outfile)
