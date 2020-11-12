# Run with --help or -h flag for help.
# Written 04/07/2020 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-N", "--infile"), type="character", default="",
		help="File assigning numbers to chromosomes (downloaded from ENA) [default= %default]", metavar="character"), 
  make_option(c("-O", "--outfile"), type="character", default="", 
		help="Enrichment file: it will be used as input, and it will be overwritten [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$infile)) {
  stop("WARNING: No infile specified with '-N' flag.")
} else {  cat ("Infile ", opt$infile, "\n")
  infile <- opt$infile  
  }

if (is.null(opt$outfile)) {
  stop("WARNING: No outfile specified with '-O' flag.")
} else {  cat ("Outfile ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }

add.chr.number<-function(infile,outfile)
{
library(data.table)
library(openxlsx)
indata<-fread(outfile,data.table=F)
if(sum(names(indata)%in%c("chr#","Chr#"))!=0) 
{
cat("All done!\n")
return()
}
numdata<-fread(infile,data.table=F)
numdata<-numdata[,c("accession","sequence-name")]
setnames(numdata,c("chr","chr#"))
indata<-merge(numdata,indata,by="chr",sort=F)
write.table(indata,outfile,sep="\t",row.names=F,quote=F)
write.xlsx(indata,gsub(".txt",".xlsx",outfile))
}
add.chr.number(infile=infile,outfile=outfile)
