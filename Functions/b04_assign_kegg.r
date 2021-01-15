# Run with --help or -h flag for help.
# Written 04/29/2019 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--infile"), type="character", default="",
		help="Formatted GO file [default= %default]", metavar="character"), 
  make_option(c("-O", "--outfile"), type="character", default="", 
		help="Formatted GO+KEGG file [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$infile)) {
  stop("WARNING: No infile specified with '-G' flag.")
} else {  cat ("Infile ", opt$infile, "\n")
  infile <- opt$infile  
  }

if (is.null(opt$outfile)) {
  stop("WARNING: No outfile specified with '-G' flag.")
} else {  cat ("Outfile ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }

assign.kegg<-function(infile,outfile)
{
library(data.table)
#library(KEGGprofile)
library(KEGGREST)
#Get all organisms present in KEGG and their code
#fullorg<-keggList("organism")
indata<-fread(infile,data.table=F)
indata$KEGG<-NA
for(aaa in 1:nrow(indata))
{
	cat("Record", aaa, "out of",nrow(indata),"\n")	
	inlist<-unlist(strsplit(indata$Gene_id[aaa],";"))
	tkegg<-rep("",length(inlist))
	for(bbb in 1:length(tkegg))
	{
	pino<-keggConv("genes", paste("uniprot",inlist[bbb],sep=":"))
	if(length(pino)<1) next
	tkegg[bbb]<-paste(pino,collapse=";")
	}
	indata$KEGG[aaa]<-paste(tkegg,collapse=";")
}
write.table(indata,outfile,sep="\t",row.names=F,quote=F)

#ciccio<-keggConv("uniprot","ath")
#Grep of the uniprot works!
#grep("Q8L6Y7",ciccio,value=T)

#This also works
#ciccio<-keggConv("ath","uniprot:Q8L6Y7")

#Superciccio is the best. We don't even need to know which one is the organism
#superciccio<-keggConv("genes", "uniprot:Q8L6Y7")

}
assign.kegg(infile=infile,outfile=outfile)
