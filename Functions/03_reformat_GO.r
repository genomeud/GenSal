# Run with --help or -h flag for help.
# Written 12/10/2019 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--infile"), type="character", default="",
				help="Original GO file [default= %default]", metavar="character"), 
  make_option(c("-O", "--out"), type="character", default="", 
				help="output file (contains reofrmatted GO/genes association) [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$infile)) {
  stop("WARNING: No (modified) DE file '-D' flag.")
} else {  cat ("DE file ", opt$infile, "\n")
  infile <- opt$infile  
  }

if (is.null(opt$out)) {
  stop("WARNING: No output file '-O' flag.")
} else {  cat ("Output file is ", opt$out, "\n")
  outfile <- opt$out  
  }

reformat.go<-function(infile,outfile)
{
library(data.table)
go<-fread(infile,data.table=F)
go_small<-go[,c("V2","V5","V6")]
cat("Removing duplicates...\n")
go_small<-go_small[!duplicated(go_small$V2),]
go_small$GO<-go_small$GO_REF<-""
setnames(go_small,"V2","Gene_id")
cat("Starting GO formatting...\n")
go_small<-go_small[,c("Gene_id","GO")]
for(aaa in 1:nrow(go_small))
{
	cat("Riga",aaa,"di",nrow(go_small),"\n")
	go_small$GO[aaa]<- paste(unique(go$V5[go$V2==go_small$Gene_id[aaa]]),sep=";",collapse=";")
}
go_small<-go_small[!duplicated(go_small$Gene_id),]
write.table(go_small,outfile,quote=F,row.names=F,sep="\t")
}
reformat.go(infile=infile,outfile=outfile)
