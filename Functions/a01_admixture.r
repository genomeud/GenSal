# Copyright:	    Gabriele Magris & Fabio Marroni 2020
# Aim:              Draw admixture plot  - Q plot + Evanno
# To add:		     
# Suggestions: 
# Fixes:  

suppressPackageStartupMessages(library("optparse"))

option_list = list(
  make_option(c("-i", "--indir"), action="store", default=NULL, type='character', 
              dest="indir", help="Input directory where file are located  [%default]"),
  make_option(c("-o", "--outdir"), action="store", default=NULL, type='character', 
              dest="outdir", help="Output directory where file are stored  [%default]"),
  make_option(c("-p", "--prefix"), action="store", default=NULL, type='character', 
              dest="file.prefix", help="File name prefix  [%default]"),
  make_option(c("--infile"), action="store", default=NULL, type='character', 
              dest="infile", help="Cross-validation summary file path  [%default]"),
  make_option(c("-k", "--k_value"), action="store", default="5.Q", type='character',
              dest="k_value",help="K value [%default]"),
  make_option(c("-P", "--popfile"), action="store", type="character",
                help="Population file path [default= %default]"), 
  make_option(c("--step"), action="store",dest="to_do",default=NULL, type='character', help="Step of the pipeline (draw_admixture,draw_cross-validation) [%default]")
			 )

opt_parser = OptionParser(option_list=option_list,description="\nDraw admixture and cross-validation plot")
opt = parse_args(OptionParser(option_list=option_list,description="\nDraw admixture and cross-validation plot"))


if(is.null(opt$k_value)) {
	print_help(opt_parser)
    stop("At least one argument must be supplied", call.=FALSE)
}

indir=opt$indir
outdir=opt$outdir
file.prefix=opt$file.prefix
infile=opt$infile
popfile=opt$popfile
k_value=opt$k_value
to_do=opt$to_do


# plot dei risultati
admixture.plot.5groups_color.no_low_cov.cov5.info50<-function(indir,outdir,file.prefix, k_value="1.Q",popfile,res=1600)
{
    library(data.table)
    # read the first two columns of the ped file in order to assign the population information
    tab<-fread(paste(indir,"/",file.prefix,".ped",sep=""),select=c(1,2),sep=" ",header=F, stringsAsFactors=F,data.table=F)
    colnames(tab)<-c("pop","sample.id")

    # read population information 
    pop_code <- read.table(popfile,stringsAsFactors=F,header=T,sep="\t")
    colnames(pop_code)<-c("sample.id","U","pop_name","origin","group_assignment")

    # merge the two datasets 
    # tab$pop<-NA
    total<-merge(tab,pop_code,by="sample.id",all.x=T,sort=F)
    # # store the dataframe with PCA coordinates 
    tab<-total[,c(6,1)]

    # translate the names 
    tab$group_assignment[tab$group_assignment=="marmorata"]<-"Marmoratus"
    tab$group_assignment[tab$group_assignment=="carpione"]<-"Carpione"
    tab$group_assignment[tab$group_assignment=="fario_atlantica"]<-"Atlantica"
    tab$group_assignment[tab$group_assignment=="fario_med_peninsula"]<-"Mediterranea Mainland"
    tab$group_assignment[tab$group_assignment=="fario_med_island"]<-"Mediterranea Island"
    samples<-tab$group_assignment
    # samples[samples=="FARIO|ATLANTICA"]<-"FARIO"


    infile=paste(outdir,"/",file.prefix,".",k_value, sep="")
    tbl=read.table(infile)
    outfile=paste(outdir,"/",file.prefix,".",k_value,"color_groups.no_low_cov.cov5.info50.jpeg", sep="")

    pop<- unique(sapply(strsplit(tab$group_assignment,"-"),"[",1))
    # pop[pop=="FARIO|ATLANTICA"]<-"FARIO"
    pop<-unique(pop)
    line_pos<-rep(0,length(pop))

    for (pp in 1:length(pop)) {
    # line_pos[pp]<-max(grep(paste("^",pop[pp],"$",sep=""), samples,fixed=T))
    line_pos[pp]<-max(which(samples==pop[pp]))
    }

    # modifico l'ordine dei campioni entro popolazione
    start<-c(1,line_pos[1:(length(line_pos)-1)]+1)
    end<-start+ diff(c(start, (nrow(tab)+1)))-1

    # append the samples ID 
    tbl$samples<-tab$sample.id


    tbl_ord<-tbl
    for (pp in 1:length(start)) {
        small<-tbl_ord[start[pp]:end[pp],]
        columns<-names(small)[order(colSums(small[-length(colnames(small))]),decreasing=TRUE)]
        tbl_ord[start[pp]:end[pp],]<-tbl_ord[start[pp]:end[pp],][order(tbl_ord[start[pp]:end[pp],columns[1]], tbl_ord[start[pp]:end[pp],columns[2]],decreasing=TRUE),]
    }

    # # # # write the new sorted Q table 
    # # # # append also the location information 
    # # # origin<-read.table("../../SNPRelate/salmo_trutta_id1492_population_assignment.txt",header=F,stringsAsFactors=F,sep="\t")
    # # # colnames(origin)<-c("samples","unlinked","population","origin")
    # # # # merge the 2 dataset + write it 
    # # # outtab<-merge(tbl_ord,origin,by=c("samples"),sort=F)
    # # # write.table(outtab,paste(file.prefix,k_value,"ordered.population.origin.txt", sep="") , col.names=T,row.names=F,sep="\t",quote=F)


    # define the line position
    line_pos<-line_pos+(line_pos*0.25)
    name_pos<-c(line_pos[1]/2, (diff(line_pos)/2+line_pos[-length(line_pos)]))

    titolo<-paste("ADMIXTURE (K=", gsub(".Q","",k_value), ")", sep="")

    res=res
    cex.text=0.5
    jpeg(outfile,width=16.9,height=8,units="cm",res=res, type="cairo")
    par(mar=c(1.2,2.5,1.4,0), mgp=c(1.6,0.5,0))
    # since I would like to specify the colors for each group (the same used for the PCA) -> based on Admixture results 
    color<-c("orchid2","green","blue","gray68","orangered")


    coord<-barplot(t(as.matrix(tbl_ord)), col=color, main=titolo, ylab="Ancestry", border=NA, cex.axis=0.8, cex.main=0.9, las=2, cex.lab=1, cex.names=0.9, space=0.25,tck=-0.02,font.lab=2)
    segments(line_pos, -0.1, x1 = line_pos, y1 = 1, lwd=1.3, xpd=NA)

    mtext(pop, side=1, line=0, at=name_pos, cex=cex.text)

    dev.off()
}


# 2. Draw cross-validation plot
plot_cv_evanno_validation<-function(infile=infile,outdir=outdir,file.prefix=file.prefix,res=1600)
{
	ylas=1
	cex.axis=0.6
	cex.point=0.8
	options(scipen=10000)
	res=res
	cex.lab=1
	cex.text=1.2
	cex.points=0.6

	t<-read.table(infile,header=T,fill=T)
	
    # draw CV plot 
	outfile=paste(outdir,"/",file.prefix,".evanno_cv_validation_error",".jpeg", sep="")
	jpeg(outfile,res=res,width=16,height=8,units="cm",type="cairo")
	par(mar=c(1.6,2.3,0.4,0.2),mgp=c(1.3,1,0))

	plot(t$K,t$CV_profile,ylim=c(0,0.6),axes=F,ylab=expression(bold(paste("Cross validation error (CV) ",sep=""))),xlab="",pch=19,cex=0.6,cex.axis=cex.axis,cex.lab=cex.lab)
	box(which="plot")
	lines(t$K,t$CV_profile,lwd=1,col="black")
	tickaxis_x<-seq(1,max(t$K),length.out=15)
	axis(1,at=tickaxis_x,las=1,cex.axis=cex.axis-.05,line=0,tck=-0.01,mgp=c(3,0.01,0),padj=-0.5)
	axis(2,at=axTicks(2),las=1,cex.axis=cex.axis,line=0,tck=-0.01,mgp=c(3,0.4,0))
	mtext("K",1,cex=cex.lab,las=1,font=2,line=0.7)

	dev.off()
	
}


if(to_do=="draw_admixture")
{
    admixture.plot.5groups_color.no_low_cov.cov5.info50(indir=indir,outdir=outdir,file.prefix=file.prefix, k_value=k_value,popfile=popfile)
} else if (to_do=="draw_cross-validation")
{
    plot_cv_evanno_validation(infile=infile,outdir=outdir,file.prefix=file.prefix,res=1600)
}