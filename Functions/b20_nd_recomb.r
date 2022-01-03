# Run with --help flag for help.
# Modified 08/01/2020 by Fabio Marroni

suppressPackageStartupMessages({
  library(optparse)
})


option_list = list(
  make_option(c("-T", "--totfile"), type="character", default="/projects/populus/ep/share/marroni/Gensal/salmo_trutta_id1492/tables/nd.cov5.info50.no_low_cov.txt", 
              help="File containing heterozygosity and HW statistics computed using vcftools [default= %default]", metavar="character"),
  make_option(c("-A", "--atlanticfile"), type="character", default="/projects/populus/ep/share/marroni/Gensal/salmo_trutta_id1492/tables/nd.cov5.info50.no_low_cov.fario_atlantica.txt", 
              help="File containing heterozygosity and HW statistics computed using vcftools [default= %default]", metavar="character"),
  make_option(c("-C", "--carpiofile"), type="character", default="/projects/populus/ep/share/marroni/Gensal/salmo_trutta_id1492/tables/nd.cov5.info50.no_low_cov.carpione.txt", 
              help="File containing heterozygosity and HW statistics computed using vcftools [default= %default]", metavar="character"),
  make_option(c("-M", "--marmofile"), type="character", default="/projects/populus/ep/share/marroni/Gensal/salmo_trutta_id1492/tables/nd.cov5.info50.no_low_cov.marmorata.txt", 
              help="File containing heterozygosity and HW statistics computed using vcftools [default= %default]", metavar="character"),
  make_option(c("-m", "--medmfile"), type="character", default="/projects/populus/ep/share/marroni/Gensal/salmo_trutta_id1492/tables/nd.cov5.info50.no_low_cov.fario_med_peninsula.txt", 
              help="File containing heterozygosity and HW statistics computed using vcftools [default= %default]", metavar="character"),
  make_option(c("-I", "--medislandfile"), type="character", default="/projects/populus/ep/share/marroni/Gensal/salmo_trutta_id1492/tables/nd.cov5.info50.no_low_cov.fario_med_island.txt", 
              help="File containing heterozygosity and HW statistics computed using vcftools [default= %default]", metavar="character"),
  make_option(c("-c", "--chrfile"), type="character", default="/projects/populus/ep/share/marroni/Gensal/genome/GCA_901001165.1_sequence_report.txt", 
              help="Minimum number of available chromosomes for the analysis [default= %default]", metavar="character"),
  make_option(c("-t", "--tol"), type="numeric", default=1000000, 
              help="Max distance between SNPs and recombination markers [default= %default]", metavar="character"),
  make_option(c("-r", "--recombfile"), type="character", default="/projects/populus/ep/share/marroni/Gensal/salmo_trutta_id1492/recomb_rate/recomb_trutta.tsv", 
              help="File containing heterozygosity and HW statistics computed using vcftools [default= %default]", metavar="character"),
  make_option(c("-O", "--outfile"), type="character", default="/projects/populus/ep/share/marroni/Gensal/salmo_trutta_id1492/recomb_rate/recomb_trutta_nd.tsv", 
              help="Output graph [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

  if (is.null(opt$totfile)) {
  stop("WARNING: No totfile specified with '-T' flag.")
} else {  cat ("totfile is", opt$totfile, "\n")
  totfile <- opt$totfile  
  }

  if (is.null(opt$atlanticfile)) {
  stop("WARNING: No atlanticfile specified with '-H' flag.")
} else {  cat ("atlanticfile is", opt$atlanticfile, "\n")
  atlanticfile <- opt$atlanticfile  
  }

  if (is.null(opt$carpiofile)) {
  stop("WARNING: No carpiofile specified with '-H' flag.")
} else {  cat ("carpiofile is", opt$carpiofile, "\n")
  carpiofile <- opt$carpiofile  
  }

  if (is.null(opt$marmofile)) {
  stop("WARNING: No marmofile specified with '-H' flag.")
} else {  cat ("marmofile is", opt$marmofile, "\n")
  marmofile <- opt$marmofile  
  }

  if (is.null(opt$medmfile)) {
  stop("WARNING: No medmfile specified with '-H' flag.")
} else {  cat ("medmfile is", opt$medmfile, "\n")
  medmfile <- opt$medmfile  
  }

  if (is.null(opt$medislandfile)) {
  stop("WARNING: No medislandfile specified with '-H' flag.")
} else {  cat ("medislandfile is", opt$medislandfile, "\n")
  medislandfile <- opt$medislandfile  
  }

  if (is.null(opt$chrfile)) {
  stop("WARNING: No chrfile specified with '-H' flag.")
} else {  cat ("chrfile is", opt$chrfile, "\n")
  chrfile <- opt$chrfile  
  }

  if (is.null(opt$recombfile)) {
  stop("WARNING: No recombfile specified with '-H' flag.")
} else {  cat ("recombfile is", opt$recombfile, "\n")
  recombfile <- opt$recombfile  
  }

  if (is.null(opt$tol)) {
  stop("WARNING: No tol specified with '-H' flag.")
} else {  cat ("tol is", opt$tol, "\n")
  tol <- opt$tol  
  }

  if (is.null(opt$outfile)) {
  stop("WARNING: No outfile specified with '-H' flag.")
} else {  cat ("outfile is", opt$outfile, "\n")
  outfile <- opt$outfile  
  }

compute_nd<-function()
{
library(data.table)
totnd<-fread(totfile,data.table=F)
carpiond<-fread(carpiofile,data.table=F)
atlnd<-fread(atlanticfile,data.table=F)
marmond<-fread(marmofile,data.table=F)
medmnd<-fread(medmfile,data.table=F)
medind<-fread(medislandfile,data.table=F)
chr<-fread(chrfile,data.table=F)
#Only keep needed info from chr conversion file 
chr<-chr[,c("accession","sequence-name")]
totnd$CHRNUM<-chr$"sequence-name"[match(totnd$CHROM,chr$accession)]
carpiond$CHRNUM<-chr$"sequence-name"[match(carpiond$CHROM,chr$accession)]
atlnd$CHRNUM<-chr$"sequence-name"[match(atlnd$CHROM,chr$accession)]
marmond$CHRNUM<-chr$"sequence-name"[match(marmond$CHROM,chr$accession)]
medmnd$CHRNUM<-chr$"sequence-name"[match(medmnd$CHROM,chr$accession)]
medind$CHRNUM<-chr$"sequence-name"[match(medind$CHROM,chr$accession)]
#Merging is no more required, I changed strategy
#allnd<-merge(farmnd,rivernd,by=c("CHRNUM","POS"),suffixes=c("_farm","_river"),sort=F)
recomb<-fread(recombfile,data.table=F)
recomb$MM_Het<-recomb$MI_Het<-recomb$MA_Het<-recomb$CA_Het<-recomb$AT_Het<-recomb$Het<-NA
#Remove high recombination rates, as Maeva did in her study
recomb<-recomb[!is.na(recomb$"recombination_rate_cM/Mb"),]
recomb<-recomb[recomb$"recombination_rate_cM/Mb"<=2,]
#Assign nd to recomb positions
for(aaa in 1:nrow(recomb))
{
stot<-totnd[totnd$CHRNUM%in%recomb$chr_num[aaa],]
stot<-stot$nd[stot$POS<=recomb$end[aaa]+tol&stot$POS>=recomb$start[aaa]-tol]
recomb$Het[aaa]<-mean(stot)
satl<-atlnd[atlnd$CHRNUM%in%recomb$chr_num[aaa],]
satl<-satl$nd[satl$POS<=recomb$end[aaa]+tol&satl$POS>=recomb$start[aaa]-tol]
recomb$AT_Het[aaa]<-mean(satl)
scarp<-carpiond[carpiond$CHRNUM%in%recomb$chr_num[aaa],]
scarp<-scarp$nd[scarp$POS<=recomb$end[aaa]+tol&scarp$POS>=recomb$start[aaa]-tol]
recomb$CA_Het[aaa]<-mean(scarp)
smarm<-marmond[marmond$CHRNUM%in%recomb$chr_num[aaa],]
smarm<-smarm$nd[smarm$POS<=recomb$end[aaa]+tol&smarm$POS>=recomb$start[aaa]-tol]
recomb$MA_Het[aaa]<-mean(smarm)
smm<-medmnd[medmnd$CHRNUM%in%recomb$chr_num[aaa],]
smm<-smm$nd[smm$POS<=recomb$end[aaa]+tol&smm$POS>=recomb$start[aaa]-tol]
recomb$MM_Het[aaa]<-mean(smm)
smi<-medind[medind$CHRNUM%in%recomb$chr_num[aaa],]
smi<-smi$nd[smi$POS<=recomb$end[aaa]+tol&smi$POS>=recomb$start[aaa]-tol]
recomb$MI_Het[aaa]<-mean(smi)
}
pdf(gsub(".tsv",".pdf",outfile))
par(mfrow=c(2,3))
par(mar=c(4.1,4.1,1.4,1))
pp<-cor.test(recomb$"recombination_rate_cM/Mb",recomb$nd)
plot(recomb$"recombination_rate_cM/Mb",recomb$nd,pch=21,ylim=c(0,0.5),ylab="Heterozygosity",xlab="",main="Total",cex=0.8,lwd=0.9,bg="white")
abline(lm(recomb$nd~recomb$"recombination_rate_cM/Mb"))
text(x=1.5,y=0.45,labels=bquote(R^2==~.(signif(pp$estimate,2))))
text(x=1.5,y=0.4,labels=bquote("P ="~.(signif(pp$p.value,3))))
pp<-cor.test(recomb$"recombination_rate_cM/Mb",recomb$AT_nd)
plot(recomb$"recombination_rate_cM/Mb",recomb$AT_nd,pch=19,ylim=c(0,0.5),ylab="",xlab="",main="AT",col="green")
abline(lm(recomb$AT_nd~recomb$"recombination_rate_cM/Mb"))
text(x=1.5,y=0.45,labels=bquote(R^2==~.(signif(pp$estimate,2))))
text(x=1.5,y=0.4,labels=bquote("P ="~.(signif(pp$p.value,3))))
pp<-cor.test(recomb$"recombination_rate_cM/Mb",recomb$CA_nd)
plot(recomb$"recombination_rate_cM/Mb",recomb$CA_nd,pch=19,ylim=c(0,0.5),ylab="",xlab="",main="CA",col="orchid2")
abline(lm(recomb$CA_nd~recomb$"recombination_rate_cM/Mb"))
text(x=1.5,y=0.45,labels=bquote(R^2==~.(signif(pp$estimate,2))))
text(x=1.5,y=0.4,labels=bquote("P ="~.(signif(pp$p.value,3))))
pp<-cor.test(recomb$"recombination_rate_cM/Mb",recomb$MA_nd)
plot(recomb$"recombination_rate_cM/Mb",recomb$MA_nd,pch=19,ylim=c(0,0.5),ylab="Heterozygosity",xlab="Recombination rate [cM/Mb]",main="MA",col="gray68")
abline(lm(recomb$MA_nd~recomb$"recombination_rate_cM/Mb"))
text(x=1.5,y=0.45,labels=bquote(R^2==~.(signif(pp$estimate,2))))
text(x=1.5,y=0.4,labels=bquote("P ="~.(signif(pp$p.value,3))))
pp<-cor.test(recomb$"recombination_rate_cM/Mb",recomb$MI_nd)
plot(recomb$"recombination_rate_cM/Mb",recomb$MI_nd,pch=19,ylim=c(0,0.5),ylab="",xlab="Recombination rate [cM/Mb]",main="MI",col="blue")
abline(lm(recomb$MI_nd~recomb$"recombination_rate_cM/Mb"))
text(x=1.5,y=0.45,labels=bquote(R^2==~.(signif(pp$estimate,2))))
text(x=1.5,y=0.4,labels=bquote("P ="~.(signif(pp$p.value,3))))
pp<-cor.test(recomb$"recombination_rate_cM/Mb",recomb$MM_nd)
plot(recomb$"recombination_rate_cM/Mb",recomb$MM_nd,pch=19,ylim=c(0,0.5),ylab="",xlab="Recombination rate [cM/Mb]",main="MM",col="orangered")
abline(lm(recomb$MM_nd~recomb$"recombination_rate_cM/Mb"))
text(x=1.5,y=0.45,labels=bquote(R^2==~.(signif(pp$estimate,2))))
text(x=1.5,y=0.4,labels=bquote("P ="~.(signif(pp$p.value,3))))
dev.off()
write.table(recomb,outfile,sep="\t",row.names=F,quote=F)
}



count.overlap<-function(start.1=1,end.1=4845859,start.2=1,end.2=84945)
{
	max.start<-max(start.1,start.2)
	min.end<-min(end.1,end.2)
	overlap<-max(0,(1+min.end-max.start))
	overlap
}

count.overlap.vector<-function(x)
{
	max.start<-max(x[1],x[3])
	min.end<-min(x[2],x[4])
	overlap<-max(0,(1+min.end-max.start))
	overlap
}



compute_nd()
