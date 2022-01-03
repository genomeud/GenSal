# Run with --help flag for help.
# Modified 08/01/2020 by Fabio Marroni
#We compared atlantica and marmorata transcripts. So the transcript "atlantica private" are those absent in marmorata.

suppressPackageStartupMessages({
  library(optparse)
})


option_list = list(
  make_option(c("-I", "--infile"), type="character", default="/projects/populus/ep/share/marroni/Gensal/salmo_trutta_id1492/recomb_rate/recomb_table.csv",
              help="Modified version of the file produced as Table S1 by Leitwein et al (the modficiation is addition of sequence names)", metavar="character"),
  make_option(c("-B", "--blastfile"), type="character", default="/projects/populus/ep/share/marroni/Gensal/salmo_trutta_id1492/recomb_rate/tag_mapped.out",
              help="Fasta output file name", metavar="character"),
  make_option(c("-G", "--graphfile"), type="character", default="/projects/populus/ep/share/marroni/Gensal/salmo_trutta_id1492/recomb_rate/recomb_trutta.pdf",
              help="Graph output file name", metavar="character"),
  make_option(c("-O", "--outfile"), type="character", default="/projects/populus/ep/share/marroni/Gensal/salmo_trutta_id1492/recomb_rate/recomb_trutta.tsv", 
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

  if (is.null(opt$blastfile)) {
  stop("WARNING: No blastfile specified with '-B' flag.")
} else {  cat ("blastfile is ", opt$blastfile, "\n")
  blastfile <- opt$blastfile  
  }

  if (is.null(opt$graphfile)) {
  stop("WARNING: No graphfile specified with '-B' flag.")
} else {  cat ("graphfile is ", opt$graphfile, "\n")
  graphfile <- opt$graphfile  
  }

  if (is.null(opt$outfile)) {
  stop("WARNING: No outfile specified with '-O' flag.")
} else {  cat ("outfile is ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }

tags_to_fasta<-function(infile,outfile,blastfile,graphfile)
{
library(data.table)
inseq<-fread(infile,data.table=F)
bres<-fread(blastfile)
#I use data table method to aggregate, thus best is a data.table
#I also did to try and start data.tables a bit
best<-bres[ , .(which.min(V10), which.max(V11)), by = .(V1)]
setnames(best,c("Name","best.e","best.bits"))
best<-best[best$best.e==best$best.bits,]
#Now best contains all the sequences for which evalue and bits are consistent in indicating the best hit
#If they are not consistent they are lost (however, for us this is not the case, so I could avoid doing this, but it was also a check)
best$start<-NA
best$end<-NA
best$chr<-NA
best$chr_long<-NA
#Stupid loop. I get the chromosome names for the unambiguous best hits
#I am sure this could have been done without a loop
for(aaa in 1:nrow(best))
{
best$chr[aaa]<-bres$V2[bres$V1==best$Name[aaa]][best$best.e[aaa]]
best$chr_long[aaa]<-bres$V12[bres$V1==best$Name[aaa]][best$best.e[aaa]]
best$start[aaa]<-bres$V8[bres$V1==best$Name[aaa]][best$best.e[aaa]]
best$end[aaa]<-bres$V9[bres$V1==best$Name[aaa]][best$best.e[aaa]]
}
#Data.table method to remove columns
best<-best[,c("best.e","best.bits"):=NULL]
#Merge the blast location on salmo trutta with the linkage map
final<-merge(inseq,best,by.x="seqnames",by.y="Name",all=T,sort=F)
final$chr_num<-trimws(unlist(lapply(strsplit(final$chr_long,":"),"[",2)))
final$chr_long<-NULL
write.table(final,outfile,row.names=F,quote=F,sep="\t")

#Plot the correspondence between LG and CHR (we will probably need a loop)
myLG<-unique(final$consensus_LG)
pdf(graphfile)
par(mfrow=c(2,2))
for(LG in 1:length(myLG))
{
small<-final[final$consensus_LG%in%myLG[LG],]
#Select the chromosomes that is more often associated with the LD. Disregard the discordant chromosomes
mytab<-table(small$consensus_LG,small$chr_num)
mynames<-colnames(mytab)
selchr<-mynames[which.max(mytab)]
small<-small[small$chr_num%in%selchr,]
plot(small$position_cM,small$start,pch=19,xlab="cM",ylab="bp",main=paste("LG",myLG[LG],", Chr",selchr))
}
dev.off()
}



tags_to_fasta(infile=infile,outfile=outfile,blastfile=blastfile,graphfile=graphfile)
