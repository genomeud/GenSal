# Run with --help flag for help.
# Modified 08/01/2020 by Fabio Marroni

suppressPackageStartupMessages({
  library(optparse)
})


option_list = list(
  make_option(c("-H", "--hwefile"), type="character", default="", 
              help="File containing heterozygosity and HW statistics computed using vcftools [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

  if (is.null(opt$hwefile)) {
  stop("WARNING: No hwefile specified with '-H' flag.")
} else {  cat ("hwefile is", opt$hwefile, "\n")
  hwefile <- opt$hwefile  
  }


compute_Fis<-function(hwefile)
{
library(data.table)
hwe<-fread(hwefile,data.table=F)
hwe$obs<-as.numeric(unlist(lapply(strsplit(hwe$"OBS(HOM1/HET/HOM2)","/"),"[",2)))
hwe$h1<-as.numeric(unlist(lapply(strsplit(hwe$"OBS(HOM1/HET/HOM2)","/"),"[",1)))
hwe$h2<-as.numeric(unlist(lapply(strsplit(hwe$"OBS(HOM1/HET/HOM2)","/"),"[",3)))
hwe$obs<-as.numeric(unlist(lapply(strsplit(hwe$"OBS(HOM1/HET/HOM2)","/"),"[",2)))
hwe$exp<-as.numeric(unlist(lapply(strsplit(hwe$"E(HOM1/HET/HOM2)","/"),"[",2)))
hwe$Fis<-1-(hwe$obs/hwe$exp)
hwe$het_obs<-hwe$obs/(hwe$h1+hwe$obs+hwe$h2)
hwe$het_exp<-hwe$exp/(hwe$h1+hwe$obs+hwe$h2)
hwe$Fis2<-1-(hwe$het_obs/hwe$het_exp)
cat("Observed heterozygosity\n")
print(summary(hwe$het_obs))
cat("Expected heterozygosity\n")
print(summary(hwe$het_exp))
cat("Fis\n")
print(summary(hwe$Fis))
cat("Fis2\n")
print(summary(hwe$Fis2))
#I just wanted to check if inbreeding measured as Fis is high as the inbreeding coefficient computed by SNPRelate.
#I saw that the order of magnitude is the same, and II am happy with it. I am not talking about this in the paper.
#I don't even write the results to file, just to screen.
}


compute_Fis(hwefile=hwefile)
