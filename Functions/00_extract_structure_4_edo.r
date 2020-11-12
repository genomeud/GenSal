# Run with --help flag for help.
# Modified 08/01/2020 by Fabio Marroni

suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-X", "--excelfile"), type="character", default="",
              help="Input excel file, basically, old version of Table S1", metavar="character"),
  make_option(c("-A", "--admixtout"), type="character", default="", 
              help="Admixture output file [default= %default]", metavar="character"),
  make_option(c("-O", "--outfile"), type="character", default="", 
              help="output file name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

  if (is.null(opt$excelfile)) {
  stop("WARNING: No excelfile specified with '-X' flag.")
} else {  cat ("excelfile is ", opt$excelfile, "\n")
  excelfile <- opt$excelfile  
  }

  if (is.null(opt$admixtout)) {
  stop("WARNING: No admixtout specified with '-A' flag.")
} else {  cat ("admixtout is ", opt$admixtout, "\n")
  admixtout <- opt$admixtout  
  }

  if (is.null(opt$outfile)) {
  stop("WARNING: No outfile specified with '-O' flag.")
} else {  cat ("outfile is ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }

attach_admixture<-function(excelfile,pedfile,admixtout,outfile)
{
library("data.table")
library(openxlsx)
metad<-read.xlsx(excelfile)
admixture<-fread(admixtout,data.table=F)
admixture<-admixture[,c("samples","V1","V2","V3","V4","V5")]
fullmetad<-merge(metad,admixture,by.x="sample_ID",by.y="samples",all=T,sort=F)
#fullmetad$uniquepop<-1-diversity(fullmetad[,c("V1","V2","V3","V4","V5")])
fullmetad$maxQ<-apply(fullmetad[,c("V1","V2","V3","V4","V5")],1,max)
setnames(fullmetad,c("V1","V2","V3","V4","V5"),c("Pop1","Pop2","Pop3","Pop4","Pop5"))
write.xlsx(fullmetad,outfile)
}

attach_admixture(excelfile=excelfile,pedfile=pedfile,admixtout=admixtout,outfile=outfile)
