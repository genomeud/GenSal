# Run with --help flag for help.
# Modified 08/01/2020 by Fabio Marroni
#We compared atlantica and marmorata transcripts. So the transcript "atlantica private" are those absent in marmorata.

suppressPackageStartupMessages({
  library(optparse)
})


option_list = list(
  make_option(c("-I", "--infile"), type="character", default="",
              help="File produced as Table S1 by Leitwein et al (https://doi.org/10.1534/g3.116.038497)", metavar="character"),
  make_option(c("-F", "--fasta"), type="character", default="",
              help="Fasta output file name", metavar="character"),
  make_option(c("-O", "--outfile"), type="character", default="", 
              help="Output file name; basically input data with sequence names added [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

  if (is.null(opt$infile)) {
  stop("WARNING: No infile specified with '-I' flag.")
} else {  cat ("infile is ", opt$infile, "\n")
  infile <- opt$infile  
  }

  if (is.null(opt$fasta)) {
  stop("WARNING: No fasta specified with '-F' flag.")
} else {  cat ("fasta is ", opt$fasta, "\n")
  fasta <- opt$fasta  
  }

  if (is.null(opt$outfile)) {
  stop("WARNING: No outfile specified with '-O' flag.")
} else {  cat ("outfile is ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }

tags_to_fasta<-function(infile,outfile,fasta)
{
library(openxlsx)
library(seqinr)
inseq<-read.xlsx(infile)
inseq$seqnames<-paste("LG",inseq$consensus_LG,"ST",inseq$chromosome,inseq$physical_position,sep="_")
write.fasta(as.list(inseq$tag),names=inseq$seqnames,as.string=F,file.out=fasta)
write.csv(inseq,outfile,row.names=F,quote=F)
}



tags_to_fasta(infile=infile,outfile=outfile,fasta=fasta)
