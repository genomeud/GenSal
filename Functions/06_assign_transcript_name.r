# Run with --help or -h flag for help.
# Written 04/07/2020 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--infile"), type="character",  default="",
		help="Formatted GO+KEGG+KEGG pathways file [default= %default]", metavar="character"), 
  make_option(c("-B", "--blastfile"), type="character",  default="",
		help="Blast tabular output file  [default= %default]", metavar="character"), 
  make_option(c("-O", "--outfile"), type="character", default="", 
		help="Formatted GO+KEGG+KEGG pathways file containing BLAST match trancript name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$infile)) {
  stop("WARNING: No infile specified with '-G' flag.")
} else {  cat ("Infile ", opt$infile, "\n")
  infile <- opt$infile  
  }

if (is.null(opt$blastfile)) {
  stop("WARNING: No blastfile specified with '-B' flag.")
} else {  cat ("blastfile ", opt$blastfile, "\n")
  blastfile <- opt$blastfile  
  }

if (is.null(opt$outfile)) {
  stop("WARNING: No outfile specified with '-G' flag.")
} else {  cat ("Outfile ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }

assign.transcript.name<-function(infile,blastfile,outfile)
{
library(data.table)
indata<-fread(infile,data.table=F)
bdata<-fread(blastfile,data.table=F)
bdata$uniprot<-unlist(lapply(strsplit(bdata$V2,"\\|"),"[",4))
bdata$uniprot<-unlist(lapply(strsplit(bdata$uniprot,"\\."),"[",1))
bdata<-bdata[,c("V1","uniprot")]
setnames(bdata,"V1","transcript")
bdata<-bdata[!duplicated(bdata),]
findata<-merge(bdata,indata,by.x="uniprot",by.y="Gene_id",sort=F)
write.table(findata,outfile,sep="\t",row.names=F,quote=F)
}
assign.transcript.name(infile=infile,blastfile=blastfile,outfile=outfile)
