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
  make_option(c("-C", "--convfile"), type="character", default="", 
              help="File for converting names of different genome releases [default= %default]", metavar="character"),
  make_option(c("-G", "--gffile"), type="character", default="", 
              help="GFF file [default= %default]", metavar="character"),
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

  if (is.null(opt$convfile)) {
  stop("WARNING: No convfile specified with '-C' flag.")
} else {  cat ("convfile is ", opt$convfile, "\n")
  convfile <- opt$convfile  
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

transcript_regions_zygosity<-function(ZHpfile,threshold,outfile,convfile,gffile)
{
library("data.table")
ZHp<-fread(ZHpfile,data.table=F)
#Remove some lines with NAs value in ZHp
ZHp<-na.omit(ZHp)
#Keep a copy of all ZHp values, we will need to compute some statistics
allZHp<-ZHp
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

#Merge the RefSeq gene names (because those are present in gff)
myconv<-fread(convfile,data.table=F)
myconv<-myconv[,c("RefSeq","INSDC")]
mysummary<-merge(mysummary,myconv,by.x="chr",by.y="INSDC",all.x=T,all.y=F,sort=F)

#Read the gff.
#I only keep the genes, we dont' care about exons
mygff<-scan(gffile,what="",sep="\n")
mygff<-mygff[substr(mygff,1,1)!="#"]
ciccio<-strsplit(mygff,"\t")
mygff<-cbind(unlist(lapply(ciccio,"[",1)),unlist(lapply(ciccio,"[",3)),unlist(lapply(ciccio,"[",4)),unlist(lapply(ciccio,"[",5)),unlist(lapply(ciccio,"[",9)))
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
	#Identify genes overlapping the region with high (or low) ZHp
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


transcript_regions_zygosity(ZHpfile=ZHpfile,threshold=threshold,outfile=outfile,convfile=convfile,gffile=gffile)
