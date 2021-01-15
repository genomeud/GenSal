# Copyright:	    Gabriele Magris 2020
# Aim:              Perform ZHp equation + draw genome-wide manhattan plot 
# To add:		     
# Suggestions: 
# Fixes:  

suppressPackageStartupMessages(library("optparse"))

option_list = list(
  make_option(c("-d", "--wdir"), action="store", default=NULL, type='character', dest="wd", help="Directory where file are located  [%default]"),
  make_option(c("-i", "--infile"), action="store", default=NULL, type='character', help="Name of input file  [%default]"),
  make_option(c("--step"), action="store",dest="to_do",default=NULL, type='character', help="Step of the pipeline (ZHp_equation,draw_plot) [%default]")
			 )
             
opt_parser = OptionParser(option_list=option_list,description="\nCalculate the ZHp equation [option ZHp_equation]\nDraw manhattan plot [option draw_plot]")
opt = parse_args(OptionParser(option_list=option_list,description="\nCalculate the ZHp equation [option ZHp_equation]\nDraw manhattan plot [option draw_plot]"))


if(is.null(opt$infile)) {
	print_help(opt_parser)
    stop("At least one argument must be supplied", call.=FALSE)
}

wd=opt$wd
infile=opt$infile
to_do=opt$to_do


zhp_equation<-function(wd="", infile="")
{
    # read the df in R and perform the Hp equation 
    library(data.table)
    setwd(wd)

    inf<-fread(infile,data.table=F)
    # loop across each rows and calculate Hp 
    inf$Hp<-NA
    out<-inf[1,]
    for(row in 1:nrow(inf))
    {
        # need to differentiate between situations in which we do not have any SNP count in the window from real counts 
        # set NA to these windows 
        if(inf$reference[row]==0 && inf$alternative[row]==-1)
        {
            # inf$Hp[row]<-NA
            # browser()
            next()
        } else {
            inf$Hp[row]<-(2*(inf$reference[row])*(inf$alternative[row]))/(inf$reference[row]+inf$alternative[row])^2
        }
        # tmp<-inf[row,]
        # tmp$Hp<-((2*tmp$reference)*(2*tmp$alternative))/(tmp$reference+tmp$alternative)^2
        # # create the output df 
        # out<-rbind(out,tmp)
    }
    # now Hp is scaled -> Z transformed 
    inf$ZHp<-scale(inf$Hp,center=T,scale=T)
    # write the output file 
    write.table(inf,gsub(".txt",".ZHp.txt",basename(infile)),quote=F,sep="\t",row.names=F)
}

plot_gw_manhattan_ZHp<-function(wd,infile)
{
    # draw now a plot of the ZHp results divided in chromosomes 
    library(data.table)
    setwd(wd)

    # set plot paramethers 
    res<-1600 
    scarto=3.5
    min.y=-5
    max.y=5
    cex.axis=0.7
    cex.lab=0.9
    cex.point=0.2
    ylas=2 
    colour=c("forestgreen","navyblue")

    # start drawing the plot 
    outfile<-paste(gsub(".txt","",basename(infile)),".jpeg",sep="")
    dir.create("../plot/",showWarnings=F,recursive=T)
    jpeg(paste("../plot/",outfile,sep=""),width=16,height=6,units="cm",res=res,type="cairo")

    par(mar=c(3,2,0.1,0.2),mgp=c(1.9,1.2,.6),oma=c(0.7,0,0.1,0),mfrow=c(1,1))


    # read input file 
    inf<-fread(infile,data.table=F)
    
    nd_red_cor<-inf

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
            nd_red_cor[nd_red_cor$index == i, ]$new_pos = nd_red_cor[nd_red_cor$index == i, ]$start/1e+06
        }
        else {
            lastbase = lastbase + tail(subset(nd_red_cor, index == 
              i - 1)$start/1e+06, 1)+scarto
            nd_red_cor[nd_red_cor$index == i, ]$new_pos = nd_red_cor[nd_red_cor$index == i, ]$start/1e+06 + 
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
    plot(nd_red_cor$new_pos,nd_red_cor$ZHp,
            xlim=c(xmin,xmax),ylim=c(min.y,max.y),xlab="",ylab="",pch=20,
            cex=0.7,col="red",xaxt = "n",yaxt = "n", bty = "o", xaxs = "i", yaxs = "i",main="",cex.axis=cex.axis,type="n",cex.lab=cex.lab,font.main=2,las=2,font.lab=2)
    # axis(2, at = axTicks(2), las=ylas,cex.axis=cex.axis,tck=-0.01,mgp=c(0.5,0.3,.6),line=0)
    axis(2, at = axTicks(2), las=ylas,cex.axis=cex.axis,tck=-0.01,mgp=c(0.5,0.3,.6),line=0)
    cex.lab=0.7
    mtext(expression(bold("ZH"["p"])),2,line=1,cex=cex.lab)
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
            tempo<-nd_red_cor[!is.na(nd_red_cor$dxy_pontica_sylvestris),]
            smoothingSpline = smooth.spline(tempo[tempo$index==i,]$new_pos,tempo[tempo$index==i,]$dxy_pontica_sylvestris,spar=smooth_spar)
            lines(smoothingSpline,col='gray45', lwd=smooth.size)
        }
        # draw points 
        points(nd_red_cor[nd_red_cor$index==i,]$new_pos,nd_red_cor[nd_red_cor$index==i,]$ZHp,col=col[icol],pch=19,cex=cex.point )

        abline(v=c(max(nd_red_cor[nd_red_cor$index==i,]$new_pos+1)),col="gray48",lwd=0.5)
        icol = icol + 1
    }
    # add the chromosome labels 
    mtext(labs,1,at=ticks,las=2,cex=par('cex')*(cex.axis-0.2),las=2,line=-0.0)
    dev.off()

}

if(to_do=="ZHp_equation")
{
    zhp_equation(wd=wd,infile=infile)
} else if (to_do=="draw_plot")
{
    plot_gw_manhattan_ZHp(wd=wd,infile=infile)
}