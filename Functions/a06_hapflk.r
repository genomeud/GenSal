# Copyright:	    Gabriele Magris & Fabio Marroni 2020
# Aim:              Draw hapflm p-value plot 
# To add:		     
# Suggestions: 
# Fixes:  

suppressPackageStartupMessages(library("optparse"))

option_list = list(
  make_option(c("-o", "--outdir"), action="store", default=NULL, type='character', 
              dest="outdir", help="Output directory where plot will be created [%default]"),
  make_option(c("-i", "--infile"), action="store", default=NULL, type='character', 
              dest="infile", help="Path to infile [%default]"),
  make_option(c("-c", "--chr_list"), action="store", default=NULL, type='character', 
              dest="chr_list", help="Ordered list of chromosomes [%default]")
              )

opt_parser = OptionParser(option_list=option_list,description="\nDraw hapflk p-value manhattan plot")
opt = parse_args(OptionParser(option_list=option_list,description="\nDraw hapflk p-value manhattan plot"))


if(is.null(opt$infile)) {
	print_help(opt_parser)
    stop("At least one argument must be supplied", call.=FALSE)
}


infile=opt$infile
outdir=opt$outdir
chr_list=opt$chr_list


# draw now a plot of the ZHp results divided in chromosomes 
draw_hapflk_pvalue_manhattan<-function(infile,outdir,chr_list)
{
library(data.table)


res<-1600 
scarto=3.5
min.y=0
max.y=3.5
cex.axis=0.7
cex.lab=0.9
cex.point=0.1
ylas=2 
colour=c("forestgreen","navyblue")

# start drawing the plot 
outfile<-paste(gsub(".txt","",basename(infile)),".jpeg",sep="")
dir.create(paste(outdir,sep=""),showWarnings=F,recursive=T)
jpeg(paste(outdir,outfile,sep=""),width=16,height=6,units="cm",res=res,type="cairo")
# png(outfile,width=10)
#par(mfrow=c(2,1))
par(mar=c(3,2,0.1,0.2),mgp=c(1.9,1.2,.6),oma=c(0.7,0,0.1,0),mfrow=c(1,1))

inf<-fread(infile,data.table=F)
inf$chro<-gsub("'","",gsub("b'","",inf$chr))
colnames(inf)[7]<-"chr"
colnames(inf)[2]<-"chro"
# read the file with the chromosome order 
chr_order<-read.table(chr_list,header=F,stringsAsFactors=F)
chr_order<-chr_order[,c(1,2)]
# select only the LRR 
head(chr_order)
chr_order<-chr_order[grep("LRR*",chr_order$V1),]
print(chr_order)
chr_order$order<-seq(1,nrow(chr_order),1)
colnames(chr_order)<-c("chr","length","order")

# merge the 2 datasets together 
out<-merge(inf,chr_order,by="chr",all.x=T,sort=F)
# order back the df 
out2<-out[order(out$order,out$pos),]
out2<-out2[!is.na(out2$order),]
nd_red_cor<-out2

# divide the chromosomes 
listchrom<-(unlist(unique(nd_red_cor$chr)))


nd_red_cor$index = NA
ind = 0
for (i in unique(nd_red_cor$chr)) {
    ind = ind + 1
    nd_red_cor[nd_red_cor$chr == i, ]$index = ind
}
# nd_red_cor<-nd_red_cor[order(nd_red_cor$index,nd_red_cor$start),]

lastbase = 0
ticks = NULL
# nd_red_cor$pos<-(nd_red_cor$start+nd_red_cor$end)/2
nd_red_cor$new_pos = 0
for (i in unique(nd_red_cor$index)) {
    if (i == 1) {
        nd_red_cor[nd_red_cor$index == i, ]$new_pos = nd_red_cor[nd_red_cor$index == i, ]$pos/1e+06
    }
    else {
        lastbase = lastbase + tail(subset(nd_red_cor, index == 
          i - 1)$pos/1e+06, 1)+scarto
        nd_red_cor[nd_red_cor$index == i, ]$new_pos = nd_red_cor[nd_red_cor$index == i, ]$pos/1e+06 + 
          lastbase
    }
    ticks = c(ticks, (min(nd_red_cor[nd_red_cor$index == i, ]$new_pos) + max(nd_red_cor[nd_red_cor$index == i, ]$new_pos))/2 + 1)
}
# xmax = ceiling(max(nd_red$new_pos) * 1.03)
# xmin = floor(max(nd_red$new_pos) * -0.03)
xmax = max(nd_red_cor$new_pos)+3.5
xmin = -3.5
xlabel = "Chromosome"
labs<-listchrom

#------------------------#
# START DRAWING THE PLOT # 
#------------------------#
plot(nd_red_cor$new_pos,-log10(nd_red_cor$pvalue),
        xlim=c(xmin,xmax),ylim=c(min.y,max.y),xlab="",ylab="",pch=20,
        cex=0.7,col="red",xaxt = "n",yaxt = "n", bty = "o", xaxs = "i", yaxs = "i",main="",cex.axis=cex.axis,type="n",cex.lab=cex.lab,font.main=2,las=2,font.lab=2)
# axis(2, at = axTicks(2), las=ylas,cex.axis=cex.axis,tck=-0.01,mgp=c(0.5,0.3,.6),line=0)
axis(2, at = axTicks(2), las=ylas,cex.axis=cex.axis,tck=-0.01,mgp=c(0.5,0.3,.6),line=0)
cex.lab=0.7
mtext(expression(bold("-log"["10"])*bolditalic(" p-value")*bold("")),2,line=1,cex=cex.lab)
mtext(expression(bold("Chromosome")),1,line=2.4,cex=cex.lab)
# # mtext(expression(bold("European")~bolditalic("sylvestris")),2,line=2.8,cex=cex.lab)
# # mtext(expression(bold(paste("vs. convarietas"))),2,line=1.8,cex=cex.lab)

# draw the threshold line 
# abline(h=tresh,lwd=1.5,col="gray68",lty=2)

col = rep(colour, max(unique(nd_red_cor$index)))
icol = 1

smoothed=F
for (i in unique(nd_red_cor$index)) {
    if(smoothed==TRUE){
        tempo<-nd_red_cor[!is.na(nd_red_cor$pvalue),]
        smoothingSpline = smooth.spline(tempo[tempo$index==i,]$new_pos,tempo[tempo$index==i,]$pvalue,spar=smooth_spar)
        lines(smoothingSpline,col='gray45', lwd=smooth.size)
    }
    # draw points 
    points(nd_red_cor[nd_red_cor$index==i,]$new_pos,-log10(nd_red_cor[nd_red_cor$index==i,]$pvalue),col=col[icol],pch=19,cex=cex.point )

    abline(v=c(max(nd_red_cor[nd_red_cor$index==i,]$new_pos+1)),col="gray48",lwd=0.5)
    icol = icol + 1
}
# add the chromosome labels 
mtext(labs,1,at=ticks,las=2,cex=par('cex')*(cex.axis-0.2),las=2,line=-0.0)
dev.off()
}

draw_hapflk_pvalue_manhattan(infile=infile,outdir=outdir,chr_list=chr_list)