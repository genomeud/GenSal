# Run with --help or -h flag for help.
# Written 04/07/2020 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--infile"), type="character", default="",
		help="Input file (ZHp, HapFLK or other similar format) [default= %default]", metavar="character"), 
  make_option(c("-K", "--keggfile"), type="character", default="",
		help="Formatted KEGG+GO file with transcript names [default= %default]", metavar="character"), 
  make_option(c("-O", "--outfile"), type="character", default="", 
		help="Enrichment output file [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$infile)) {
  stop("WARNING: No infile specified with '-G' flag.")
} else {  cat ("Infile ", opt$infile, "\n")
  infile <- opt$infile  
  }

if (is.null(opt$keggfile)) {
  stop("WARNING: No keggfile specified with '-K' flag.")
} else {  cat ("keggfile ", opt$keggfile, "\n")
  keggfile <- opt$keggfile  
  }

if (is.null(opt$outfile)) {
  stop("WARNING: No outfile specified with '-G' flag.")
} else {  cat ("Outfile ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }

kegg.enrich<-function(infile,keggfile,outfile,minpos=3,p.value=0.05)
{
library(data.table)
library(openxlsx)
indata<-fread(infile,data.table=F)
bdata<-fread(keggfile,data.table=F)
#Remove lines without KEGG pathway
bdata<-bdata[bdata$pKEGG_class!=""|bdata$pKEGG_map!="",]
#Create the universe of KEGG classes and maps
allkeggclass<-unlist(strsplit(bdata$pKEGG_class,"\\|"))
allkeggclass<-allkeggclass[allkeggclass!=""]
allkeggmap<-unlist(strsplit(bdata$pKEGG_map,"\\|"))
allkeggmap<-allkeggmap[allkeggmap!=""]
indata$mappval<-indata$mapbg<-indata$mapratio<-indata$keggmap<-indata$classpval<-indata$classbg<-indata$classratio<-indata$keggclass<-""
#Loop over significant results
for(aaa in 1:nrow(indata))
{
	cat("aaa=",aaa,"\n")
	genewin<-data.frame(genewin=unlist(strsplit(indata$genbank[aaa],";")))
	ncbiwin<-data.frame(genewin=unlist(strsplit(indata$gene_name[aaa],";")))
	pp<-merge(genewin,bdata,by.x="genewin",by.y="transcript")
	mykeggclass<-unlist(strsplit(pp$pKEGG_class,"\\|"))
	mykeggclass<-mykeggclass[mykeggclass!=""]
	mykeggmap<-unlist(strsplit(pp$pKEGG_map,"\\|"))
	mykeggmap<-mykeggmap[mykeggmap!=""]
	myuniqclass<-unique(mykeggclass)
	keep<-classratio<-classbg<-classp<-rep(0,length(myuniqclass))
	#if(aaa==3) browser()
	if(length(myuniqclass)>0)
	{
		for(bbb in 1:length(myuniqclass))
		{
			cat("first bbb=",bbb,"\n")
			set1<-sum(mykeggclass==myuniqclass[bbb])
			set0<-length(mykeggclass)
			bg1<-sum(allkeggclass==myuniqclass[bbb])
			bg0<-length(allkeggclass)
			classratio[bbb]<-paste(set1,set0,sep="/")
			classbg[bbb]<-paste(bg1,bg0,sep="/")
			classp[bbb]<-fisher.test(matrix(c(set1,set0,bg1,bg0),nrow=2))$p.value
			#Only write entries with at least minpos terms in the same KEGG category and with a p.value for enrichment lower than p.value 
			if(set1>=minpos & classp[bbb]<=p.value) keep[bbb]<-1
		}
		indata$keggclass[aaa]<-paste(myuniqclass[keep>0],collapse="|")
		indata$classratio[aaa]<-paste(classratio[keep>0],collapse="|")
		indata$classbg[aaa]<-paste(classbg[keep>0],collapse="|")
		indata$classpval[aaa]<-paste(classp[keep>0],collapse="|")
	}
	myuniqmap<-unique(mykeggmap)
	keep<-mapratio<-mapbg<-mapp<-rep(0,length(myuniqmap))
	if(length(myuniqmap)>0)
	{
		for(bbb in 1:length(myuniqmap))
		{
			set1<-sum(mykeggmap==myuniqmap[bbb])
			set0<-length(mykeggmap)
			bg1<-sum(allkeggmap==myuniqmap[bbb])
			bg0<-length(allkeggmap)
			mapratio[bbb]<-paste(set1,set0,sep="/")
			mapbg[bbb]<-paste(bg1,bg0,sep="/")
			mapp[bbb]<-fisher.test(matrix(c(set1,set0,bg1,bg0),nrow=2))$p.value
			#Only write entries with at least minpos terms in the same KEGG category and with a p.value for enrichment lower than p.value 
			if(set1>=minpos & mapp[bbb]<=p.value) keep[bbb]<-1

		}
		indata$keggmap[aaa]<-paste(myuniqmap[keep>0],collapse="|")
		indata$mapratio[aaa]<-paste(mapratio[keep>0],collapse="|")
		indata$mapbg[aaa]<-paste(mapbg[keep>0],collapse="|")
		indata$mappval[aaa]<-paste(mapp[keep>0],collapse="|")
	}
}
write.table(indata,outfile,sep="\t",row.names=F,quote=F)
write.xlsx(indata,gsub(".txt",".xlsx",outfile))
}
kegg.enrich(infile=infile,keggfile=keggfile,outfile=outfile)
