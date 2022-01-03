# Run with --help flag for help.
# Modified 08/01/2020 by Fabio Marroni

suppressPackageStartupMessages({
  library(optparse)
})


option_list = list(
  make_option(c("-F", "--frqfile"), type="character", default="", 
              help="File containing heterozygosity and HW statistics computed using vcftools [default= %default]", metavar="character"),
  make_option(c("-c", "--minchr"), type="numeric", default=20, 
              help="Minimum number of available chromosomes for the analysis [default= %default]", metavar="character"),
  make_option(c("-O", "--outfile"), type="character", default="", 
              help="Output file with nd for each SNP [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

  if (is.null(opt$frqfile)) {
  stop("WARNING: No frqfile specified with '-H' flag.")
} else {  cat ("frqfile is", opt$frqfile, "\n")
  frqfile <- opt$frqfile  
  }

  if (is.null(opt$minchr)) {
  stop("WARNING: No minchr specified with '-H' flag.")
} else {  cat ("minchr is", opt$minchr, "\n")
  minchr <- opt$minchr  
  }

  if (is.null(opt$outfile)) {
  stop("WARNING: No outfile specified with '-H' flag.")
} else {  cat ("outfile is", opt$outfile, "\n")
  outfile <- opt$outfile  
  }

compute_nd<-function()
{
library(data.table)
hwe<-fread(frqfile,data.table=F)
#Stupid trick to manage horrible name convention in frq files
names(hwe)<-c("CHROM","POS","N_ALLELES","N_CHR","REF","ALT")
hwe$N_ALLELES<-NULL
#Conver AF to numeric
hwe$REF<-as.numeric(unlist(lapply(strsplit(hwe$REF,":"),"[",2)))
hwe$ALT<-as.numeric(unlist(lapply(strsplit(hwe$ALT,":"),"[",2)))
#Compute nucleotide diversity as 2*p*q
hwe$nd<-2*hwe$REF*hwe$ALT
#Filter for number of chromosomes
hwe<-hwe[hwe$N_CHR>=minchr,]
write.table(hwe,outfile,sep="\t",row.names=F,quote=F)
}


compute_nd()
