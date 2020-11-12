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
  make_option(c("-C", "--convfile"), type="character", default="", 
              help="File for converting names of different genome releases [default= %default]", metavar="character"),
  make_option(c("-D", "--maxdist"), type="numeric", default=1000000, 
              help="Maximum distance for two signal to be clustered [default= %default]", metavar="character"),
  make_option(c("-S", "--minsnp"), type="numeric", default=5, 
              help="Minimum number of consecutive significant SNPs to be included in results [default= %default]", metavar="character"),
  make_option(c("-G", "--gffile"), type="character", default="", 
              help="GFF file [default= %default]", metavar="character"),
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

  if (is.null(opt$convfile)) {
  stop("WARNING: No convfile specified with '-C' flag.")
} else {  cat ("convfile is ", opt$convfile, "\n")
  convfile <- opt$convfile  
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

  if (is.null(opt$gffile)) {
  stop("WARNING: No gffile specified with '-G' flag.")
} else {  cat ("gffile is ", opt$gffile, "\n")
  gffile <- opt$gffile  
  }

  if (is.null(opt$outfile)) {
  stop("WARNING: No outfile specified with '-O' flag.")
} else {  cat ("outfile is ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }

extract_HFLK_sig<-function(infile,threshold,outfile,convfile,gffile,pvalue,maxdist,minsnp)
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

for(aaa in 2:nrow(hFLK))
{
cat("Riga",aaa,"di",nrow(hFLK),"\n")
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

#Merge the RefSeq chromosome names (because those are present in gff)
myconv<-fread(convfile,data.table=F)
myconv<-myconv[,c("RefSeq","INSDC")]
mysummary<-merge(finalFLK,myconv,by.x="chr",by.y="INSDC",all.x=T,all.y=F,sort=F)

#Read the gff using scan, because the irregular format could create problems to fread.
#Why couldn't they create a GFF with a fixed number of columns? I don't know
mygff<-scan(gffile,what="",sep="\n")
#Skip header
mygff<-mygff[substr(mygff,1,1)!="#"]
#Split by tab and only get five fixed columns
ciccio<-strsplit(mygff,"\t")
mygff<-cbind(unlist(lapply(ciccio,"[",1)),unlist(lapply(ciccio,"[",3)),unlist(lapply(ciccio,"[",4)),unlist(lapply(ciccio,"[",5)),unlist(lapply(ciccio,"[",9)))
#I only keep the genes, we dont' care about exons
newgff<-mygff[mygff[,2]=="mRNA",]
newgff<-data.frame(newgff,stringsAsFactors=F)
newgff$X2<-NULL
#Extract gene names and gene ID from the gff
newgff$genbank<-unlist(lapply(strsplit(newgff$X5,";"),"[",3))
newgff$genbank<-gsub("Genbank:","",unlist(lapply(strsplit(newgff$genbank,","),"[",2)))
#newgff$gene_name<-gsub("gene=","",unlist(lapply(strsplit(newgff$X5,";"),"[",6)))
newgff$gene_name<-unlist(lapply(strsplit(newgff$X5,"Dbxref=GeneID:"),"[",2))
newgff$gene_name<-unlist(lapply(strsplit(newgff$gene_name,","),"[",1))
newgff$description<-unlist(lapply(strsplit(newgff$X5,"product="),"[",2))
newgff$description<-unlist(lapply(strsplit(newgff$description,";"),"[",1))
newgff$description<-unlist(lapply(strsplit(newgff$description,"%"),"[",1))
#newgff$gene_ID<-gsub("Dbxref=GeneID:","",unlist(lapply(strsplit(newgff$X5,";"),"[",2)))
newgff$X5<-NULL
setnames(newgff,c("X1","X3","X4"),c("chr","start","end"))
newgff$start<-as.numeric(as.character(newgff$start))
newgff$end<-as.numeric(as.character(newgff$end))
mysummary$description<-mysummary$gene_name<-mysummary$genbank<-""
for(bbb in 1:nrow(mysummary))
{
	cat("Riga",bbb,"di",nrow(mysummary),"\n")
	sgff<-newgff[newgff$chr%in%mysummary$RefSeq[bbb],]
	sgff$ov<-NA
	#Identify genes overlapping the region with significant HFLK signal
	for(ccc in 1:nrow(sgff))
	{
	sgff$ov[ccc]<-count.overlap(sgff$start[ccc],sgff$end[ccc],mysummary$start[bbb],mysummary$end[bbb])
	}
	#Print all the genes names located in each significant window
	mysummary$genbank[bbb]<-paste(unique(sgff$genbank[sgff$ov>0]),sep=";",collapse=";")
	mysummary$gene_name[bbb]<-paste(unique(sgff$gene_name[sgff$ov>0]),sep=";",collapse=";")
	mysummary$description[bbb]<-paste(unique(sgff$description[sgff$ov>0]),sep=";",collapse=";")
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


extract_HFLK_sig(infile=infile,pvalue=pvalue,maxdist=maxdist,minsnp=minsnp,outfile=outfile,convfile=convfile,gffile=gffile)
