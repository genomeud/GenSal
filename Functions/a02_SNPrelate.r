# Run with --help or -h flag for help.
# Written 04/07/2020 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-V", "--vcffile"), type="character",
  default="salmo_trutta_id1492/stacks/populations.snps.cov5.info50.no_low_cov.vcf",help="stacks vcf file [default= %default]", metavar="character"), 
  make_option(c("-P", "--popfile"), type="character",
  default="salmo_trutta_id1492/salmo_trutta_id1492_population_assignment.5groups.txt",help="Population file path [default= %default]", metavar="character"), 
  make_option(c("-G", "--gdsfile"), type="character",
  default="salmo_trutta_id1492/SNPRelate/salmo_trutta_id1492.cov5.info50.gds",help="SNPRelate gds file [default= %default]", metavar="character"), 
  make_option(c("-B", "--between"), type="character",
  default="salmo_trutta_id1492/SNPRelate/IBD_between.pdf",help="Between population IBD graph file [default= %default]", metavar="character"), 
  make_option(c("-W", "--within"), type="character",
  default="salmo_trutta_id1492/SNPRelate/IBD_within.pdf",help="Within population IBD graph file [default= %default]", metavar="character"), 
  make_option(c("-M", "--maf"), type="numeric",
  default=0.05,help="Minimum MAF to include a SNP in analysis [default= %default]", metavar="character"), 
  make_option(c("-S", "--missingness"), type="numeric",
  default=0.05,help="Maximum missingness to include a SNP in analysisi [default= %default]", metavar="character"), 
  make_option(c("-L", "--LD"), type="numeric",
  default=0.2,help="Threshold for LD pruning [default= %default]", metavar="character"), 
  make_option(c("-I", "--ibdmatfile"), type="character",
  default="salmo_trutta_id1492/SNPRelate/salmo_trutta_id1492.cov5.info50.ibd", help="IBD matrix output file [default= %default]", metavar="character"),
  make_option(c("--step"), action="store",dest="to_do",default=NULL, type='character', help="Step of the pipeline (pca_analysis,ibd_analysis) [%default]")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")


if (is.null(opt$vcffile)) {
  stop("WARNING: No vcf file specified with '-V' flag.")
} else {  cat ("vcf file ", opt$gdsfile, "\n")
  gdsfile <- opt$gdsfile  
  }

  if (is.null(opt$popfile)) {
  stop("WARNING: No popfile specified with '-P' flag.")
} else {  cat ("popfile is", opt$popfile, "\n")
  popfile <- opt$popfile  
  }

if (is.null(opt$gdsfile)) {
  stop("WARNING: No gds file specified with '-G' flag.")
} else {  cat ("gds file ", opt$gdsfile, "\n")
  gdsfile <- opt$gdsfile  
  }
  
  if (is.null(opt$between)) {
  stop("WARNING: No between file specified with '-B' flag.")
} else {  cat ("between is", opt$between, "\n")
  between <- opt$between  
  }

  if (is.null(opt$within)) {
  stop("WARNING: No within file specified with '-W' flag.")
} else {  cat ("within is", opt$within, "\n")
  within <- opt$within  
  }

  if (is.null(opt$highcovdir)) {
  stop("WARNING: No highcov dir specified with '-H' flag.")
} else {  cat ("highcovdir is", opt$highcovdir, "\n")
  highcovdir <- opt$highcovdir  
  }

  if (is.null(opt$maf)) {
  stop("WARNING: No maf specified with '-M' flag.")
} else {  cat ("maf is", opt$maf, "\n")
  maf <- opt$maf  
  }

  if (is.null(opt$missingness)) {
  stop("WARNING: No missingness specified with '-S' flag.")
} else {  cat ("missingness is", opt$missingness, "\n")
  missingness <- opt$missingness  
  }

  if (is.null(opt$LD)) {
  stop("WARNING: No LD specified with '-L' flag.")
} else {  cat ("LD is", opt$LD, "\n")
  LD <- opt$LD  
  }

  if (is.null(opt$ibdmatfile)) {
  stop("WARNING: No ibdmatfile specified with '-I' flag.")
} else {  cat ("ibdmatfile is", opt$ibdmatfile, "\n")
  ibdmatfile <- opt$ibdmatfile  
  }
  
  if (is.null(opt$to_do)) {
  stop("WARNING: No step specified with '--step' flag.")
} else {  cat ("Step to run is", opt$to_do, "\n")
  to_do <- opt$to_do  
  }

dir.create(dirname(ibdmatfile),showWarnings=F,recursive=T)

run_pca<-function(vcffile,popfile)
{
    library(gdsfmt)
    library(SNPRelate)
	threads=4
	
	dir.create("SNPRelate",recursive=T,showWarnings=F)
	outpath<-"SNPRelate/"
	outsuffix<-"salmo_trutta_id1492"
    #######
    # PCA #
    #######
    vcf.fn<-vcffile
    snpgdsVCF2GDS(vcf.fn, "SNPRelate/salmo_trutta_id1492.cov5.info50.no_low_cov.gds", method="biallelic.only")

    # it is possible to read an already created gds object:
	(genofile <- snpgdsOpen("SNPRelate/salmo_trutta_id1492.cov5.info50.no_low_cov.gds"))	

	sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
    # exclude sylvestris + outgroups 
    # # to_exclude<-strsplit(to_remove, ",")[[1]]
    # # print(to_exclude)
    # # sample.id<-sample.id[! sample.id %in% to_exclude]


    # calculate PCA 
    pca<-snpgdsPCA(genofile,sample.id=sample.id, num.thread=threads,autosome.only=FALSE)
    # create output dataframe with the eigenvalues 
    tab <- data.frame(sample.id=pca$sample.id,EV1 = pca$eigenvect[,1], EV2=pca$eigenvect[,2],EV3 = pca$eigenvect[,3], EV4=pca$eigenvect[,4],stringsAsFactors=F)

    # read population information 
    pop_code <- read.table(popfile,stringsAsFactors=F,header=F,sep="\t")
    colnames(pop_code)<-c("sample.id","U","pop_name","origin","group_assignment")
    # merge the two datasets 
    # tab$pop<-NA
    total<-merge(tab,pop_code,by="sample.id",all.x=T,sort=F)
    # # store the dataframe with PCA coordinates 
    write.table(total,paste(outpath,"PCA_",res,".all_snp.",outsuffix,".txt",sep=""),sep="\t",quote=F,row.names=F,col.names=T)
    # total<-read.table(paste(outpath,"PCA_",res,".all_snp.",outsuffix,".",group,".txt",sep=""),header=T,stringsAsFactors=F,check.names=F)


    ########
    # PLOT #
    ########
    cex.point=0.9
    cex.axis=0.5
    cex.cross=0.7
    cex.sylvestrys=0.5
    cex.legend=0.7
    cex.text=1.8

    out_color=paste(outpath,"PCA_",res,".all_snp.",outsuffix,".final_group_assignment.removed_low_cov.cov5.info50.jpeg",sep="")
    jpeg(out_color,width=16,height=8,units="cm",res=res,type="cairo")
    par(mgp=c(1.25,0.3,0),oma=c(0.2,0.2,0.05,0.2),mar=c(2,3,1,1))

    #################
    # ASSIGN COLORS #
    #################
    col1="blue"
    col2="green"
    col3="orangered"
    col4="gray68"
    col5="orchid2"

    total$color_group<-"black"
    total$color_group[total$group_assignment=="fario_med_island"]<-col1
    total$color_group[total$group_assignment=="fario_atlantica"]<-col2
    total$color_group[total$group_assignment=="fario_med_peninsula"]<-col3
    total$color_group[total$group_assignment=="marmorata"]<-col4
    total$color_group[total$group_assignment=="carpione"]<-col5

    
    # test<-total
    color_list<-c(col1,col2,col3,col4,col5,col6,col7,col8,col9,col10)

    test<-total
 	########
	# AXES #
	########
	# set axes -> round up to nearest --> pay attention to sign
	if (min(test$EV1)<0){
		x_le<-floor(min(test$EV1)*100)/100
	} else {
		x_le<-ceiling(min(test$EV1)*100)/100
	}
	if (max(test$EV1)<0){
		x_ri<-floor(max(test$EV1)*100)/100
	} else {
		x_ri<-ceiling(max(test$EV1)*100)/100
	}
	if (max(test$EV2)<0){
		y_up<-floor(max(test$EV2)*100)/100
	} else {
		y_up<-ceiling(max(test$EV2)*100)/100
	}
	if (min(test$EV2)<0){
		y_low<-floor(min(test$EV2)*100)/100
	} else {
		y_low<-ceiling(min(test$EV2)*100)/100
	}


	############
	# VARIANCE #
	############
	# calculate % variance of each eigenvector 
    # variance proportion (%)
    pc.percent <- pca$varprop*100
    head(round(pc.percent, 2))

	pc1=1
	pc2=2
	plot(test$EV1,test$EV2,col=test$color_group,cex=cex.cross,pch=c(0),xlim=c(x_le,x_ri),ylim=c(y_low,y_up),xlab="",ylab="",font.lab=2,las=1,cex.axis=cex.axis,xaxt="n",yaxt="n",tck=-0.01)



	########
	# AXES #
	########
	#draw the x axis , since tick labels are to far from axis
	axis(1,at=axTicks(1),labels=round((axTicks(1)),2),cex.axis=cex.axis+.2,lwd=1,line=0,tck=-0.01,mgp=c(1.6,0.01,0),las=0,font.lab=2)
	# add the x axis of PC1 for Myles as text
	mtext(paste("PC",pc1," (",round(pc.percent[1],2), "%)",sep=""),1,font=2,cex=1,line=1.2)
	
	#draw the y axis , since tick labels are to far from axis
	axis(2,at=round(axTicks(2),2),,labels=round((axTicks(2)),2),cex.axis=cex.axis+.2,lwd=1,line=0,tck=-0.01,mgp=c(1.6,0.3,0),las=0,font.lab=2,las=1)
	# add the y axis of PC1 for Myles as text
	mtext(paste("PC",pc2," (",round(pc.percent[2],2), "%)",sep=""),2,font=2,cex=1,line=1.8)

	##########
	# LEGEND #
	##########
    par(family="Calibri")
    legend("topleft",legend=c(expression(paste(italic("Mediterranea Island"))),expression(paste(italic("Atlantica"))),expression(paste(italic("Mediterranea Mainland"))),expression(paste(italic("Marmoratus"))),expression(paste(italic("Carpione")))),text.col=c(c(color_list)[1:5]),ncol=1,x.intersp=0.7,text.font=c(4),cex=cex.legend,pch=c(0),col=c(c(color_list)[1:5]),bty = "n",bg =NA)

   
    dev.off()

    snpgdsClose(genofile)
  
  
}





run_ibd<-function(gdsfile,coefficient,popfile,between,within,ibdmatfile,maf,missingness,LD)
{
    library(data.table)
    library(RColorBrewer)
    library(gdsfmt)
    library(SNPRelate)
    library(dplyr)
    library(multcompView)


    #Load gds file (has to be previously saved)
    genofile <- snpgdsOpen(gdsfile)


    #Prune for LD (otherwise IBD estimates might be inflated and running time longer)
    snpset <- snpgdsLDpruning(genofile, ld.threshold=LD,autosome.only=F)
    snpset.id<-unlist(snpset)

    sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
    #Compute IBD and write IBD matrix just in case...
    ibd <- snpgdsIBDMLE(genofile, maf=maf, missing.rate=missingness, num.thread=8,snp.id=snpset.id,autosome.only=F)
    ibd.coeff <- snpgdsIBDSelection(ibd)
    ibd.coeff<-ibd.coeff[order(ibd.coeff$kinship,decreasing=T),]

    write.table(ibd.coeff,ibdmatfile,row.names=F,quote=F,sep="\t")

    #Assign populations to individuals
    mypop<-fread(popfile,data.table=F)
    mypop$final_population[mypop$final_population=="fario_atlantica"]<-"AT"
    mypop$final_population[mypop$final_population=="fario_med_island"]<-"MI"
    mypop$final_population[mypop$final_population=="fario_med_peninsula"]<-"MM"
    mypop$final_population[mypop$final_population=="carpione"]<-"CA"
    mypop$final_population[mypop$final_population=="marmorata"]<-"MA"

    ibd.coeff$POP2<-ibd.coeff$POP1<-NULL
    ibd.coeff$POP1<-mypop$final_population[match(ibd.coeff$ID1,mypop$sample_ID)]
    ibd.coeff$POP2<-mypop$final_population[match(ibd.coeff$ID2,mypop$sample_ID)]

    #Only select within population comparisons
    intra<-ibd.coeff[ibd.coeff$POP1==ibd.coeff$POP2,]

    # create also vs all population
    intra2<-ibd.coeff
    intra2$POP1<-"all"

    intra<-rbind(intra,intra2)
    #Plot within pop relatedness
    intra$POP1<-factor(intra$POP1,levels=c("MI","AT","MM","MA","CA","all"))
    #intra$ind1.id<-factor(intra$ind1.id,levels=c("Mediterranea\nIsland","Atlantica","Mediterranea\nMainland","Marmoratus","Carpione"))
    color<-c("blue","green","orangered","gray68","orchid2","white")

    # ---------------------- #
    # inbreeding coefficient #
    # ---------------------- #
    # calculate the inbdreeding coefficient
    #Prune for LD (otherwise IBD estimates might be inflated and running time longer)
    snpset <- snpgdsLDpruning(genofile, ld.threshold=LD,autosome.only=F)
    snpset.id<-unlist(snpset)

    sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

    indinb.coeff.mom.weir <- snpgdsIndInb(genofile, maf=maf, missing.rate=missingness, snp.id=snpset.id,autosome.only=F,method="mom.weir")
    indinb.coeff.mom.visscher <- snpgdsIndInb(genofile, maf=maf, missing.rate=missingness,snp.id=snpset.id,autosome.only=F,method="mom.visscher")
    indinb.coeff.mle <- snpgdsIndInb(genofile, maf=maf, missing.rate=missingness, snp.id=snpset.id,autosome.only=F,method="mle")

    # assign now the value to each individual / population <- i don't understand the order of samples -
    # create a df with sample id + inbred + population
    out<-data.frame(cbind(sample.id,indinb.coeff.mle$inbreeding),stringsAsFactors=F)
    colnames(out)[2]<-"indinb.coeff.mle"
    out$indinb.coeff.mle<-as.numeric(out$indinb.coeff.mle)
    out$POP1<-mypop$final_population[match(out$sample.id,mypop$sample_ID)]
    #Plot within pop relatedness
    out$POP1<-factor(out$POP1,levels=c("MI","AT","MM","MA","CA"))
    # create a new distribution which include all samples 
    out2=out
    out2$POP1="all"
    out<-rbind(out,out2)


    ################
    # PLOT drawing # 
    ################

    # draw a merged plot that integrate kinship + inbreeding
    res=1200
    cex.text=0.6
    cex_text=1.6
    cex.axis=0.75
    cex.significance=1
    cex.lab=1
    ylas=1
    pos=1

    outfile<-paste(dirname(within),"/kinship_within.inbreeding_coefficient.v2.jpeg",sep="")
    jpeg(outfile,width=16,height=12,units="cm",res=res, type="cairo")
    par(mar=c(0.8,3,0.2,0.1), mgp=c(1.6,0.5,0),tck=-0.03,oma=c(1.5,0.1,0.3,0.1),mfrow=c(2,1))

    ###################
    # PLOT A # whitin #
    ###################

    posi<-boxplot(intra$kinship~intra$POP1,col=color,cex.axis=0.8,ylab="",xlab="",font=2,yaxt="n",xaxt="n",ylim=c(-0.02,0.5))
    axis(2,line=0,las=ylas,cex.axis=cex.axis,mgp=c(0,0.5,0))
    mtext("Kinship coefficient",2,font=2,line=1.8,cex=cex.lab)
    # add the panel letter
    mtext(LETTERS[pos],2,line=2,at=0,padj=-8.5,adj=1,cex=cex_text,las=2,font=2)
    pos=pos+1


    # add the pairwise wilcoxon test 
    # create a dataframe with the values from the matrix 
    pp<-pairwise.wilcox.test(intra$kinship,intra$POP1,p.adjust.method="none", paired = FALSE)
    mymat<-tri_to_squ(pp$p.value)
    myletters<-multcompLetters(mymat,compare="<=",threshold=0.05,Letters=letters)
    print(myletters)
    print(posi$stats)

    # # extract the order (median) of groups
    ordine<-posi$names[order(posi$stats[3,])]
    lista<-as.list(strsplit(myletters$Letters,split=""))
    # create a new order based on the median average (obtained from boxplot data)
    intra$ordinato<-0
    for(tt in 1:length(ordine)){
        intra$ordinato[intra$POP1==ordine[tt]]<-tt
    }
    # now need to reassign the corresponding data to the names 
    pp1<-pairwise.wilcox.test(intra$kinship,intra$ordinato,p.adjust.method="none", paired = FALSE)
    # the colum and row have now the ordinato information 
    colnames(pp1$p.value)<-ordine[-length(ordine)]
    rownames(pp1$p.value)<-ordine[-1]
    mymat<-tri_to_squ(pp1$p.value)
    myletters[1]$Letters<-multcompLetters(mymat,compare="<=",threshold=0.05,Letters=letters)$Letters
    # order back as in boxplot in order to add the text to each box
    # myletters<-myletters[1]$Letters[order(names(myletters[1]$Letters))]
    myletters<-myletters[1]$Letters
    # myletters<-myletters[order(match(as.character(unique(comp1$V2)),names(myletters)))]
    myletters<-myletters[order(match(names(myletters),as.character(levels(intra$POP1))))]
    print(myletters)
    # for the y positions use 0
    # stat$stats[5,]+(letter_p_shift*(stat$stats[5,])/100)
    text(seq(1,length(posi$names),1),par('usr')[3]-((55*par('usr')[3])/100),myletters,cex=cex.significance,font=2,adj=c(0.5,0.5) )


    #############################
    # PLOT B # inbreeding (MLE) #
    #############################
    posi<-boxplot(out$indinb.coeff.mle~out$POP1,col=color,cex.axis=0.8,ylab="",xlab="",lab.font=2,yaxt="n",xaxt="n",ylim=c(0,1))
    axis(2,line=0,las=ylas,cex.axis=cex.axis,mgp=c(0,0.5,0))
    mtext("Inbreeding coefficient",2,font=2,line=1.8,cex=cex.lab)
    axis(1,at=c(1:6),tick=T,labels=FALSE,adj=0.5,padj=0.5,font=c(4),cex.axis=cex.axis)
    text(c(1:5),-0.18,c("Mediterranea\nIsland","Atlantica","Mediterranea\nMainland","Marmoratus","Carpione"),adj=0.5,padj=0.5,font=4,cex=cex.axis,xpd=NA,las=1)
    text(c(6),-0.18,c("All"),adj=0.5,padj=0.5,font=2,cex=cex.axis,xpd=NA,las=1)
    # axis(1,at=c(1:6),labels=c("Mediterranea\nIsland","Atlantica","Mediterranea\nMainland","Marmoratus","Carpione"),adj=0.5,padj=0.5,font=c(4),cex.axis=cex.axis)
    # add the panel letter
    mtext(LETTERS[pos],2,line=2,at=0,padj=-9.3,adj=1,cex=cex_text,las=2,font=2)



    # add the pairwise wilcoxon test 
    # create a dataframe with the values from the matrix 
    pp<-pairwise.wilcox.test(out$indinb.coeff.mle,out$POP1,p.adjust.method="none", paired = FALSE)
    mymat<-tri_to_squ(pp$p.value)
    myletters<-multcompLetters(mymat,compare="<=",threshold=0.05,Letters=letters)
    print(myletters)
    print(posi$stats)

    # # extract the order (median) of groups
    ordine<-posi$names[order(posi$stats[3,])]
    lista<-as.list(strsplit(myletters$Letters,split=""))
    # create a new order based on the median average (obtained from boxplot data)
    out$ordinato<-0
    for(tt in 1:length(ordine)){
        out$ordinato[out$POP1==ordine[tt]]<-tt
    }
    # now need to reassign the corresponding data to the names 
    pp1<-pairwise.wilcox.test(out$indinb.coeff.mle,out$ordinato,p.adjust.method="none", paired = FALSE)
    # the colum and row have now the ordinato information 
    colnames(pp1$p.value)<-ordine[-length(ordine)]
    rownames(pp1$p.value)<-ordine[-1]
    mymat<-tri_to_squ(pp1$p.value)
    myletters[1]$Letters<-multcompLetters(mymat,compare="<=",threshold=0.05,Letters=letters)$Letters
    # order back as in boxplot in order to add the text to each box
    # myletters<-myletters[1]$Letters[order(names(myletters[1]$Letters))]
    myletters<-myletters[1]$Letters
    # myletters<-myletters[order(match(as.character(unique(comp1$V2)),names(myletters)))]
    myletters<-myletters[order(match(names(myletters),as.character(levels(out$POP1))))]
    # for the y positions use 0
    # stat$stats[5,]+(letter_p_shift*(stat$stats[5,])/100)
    print(par('usr')[3])
    text(seq(1,length(posi$names),1),par('usr')[3]-((125*par('usr')[3])/100),myletters,cex=cex.significance,font=2,adj=c(0.5,0.5) )

    dev.off()

}



tri_to_squ<-function(x)
{
    rn<-row.names(x)
    cn<-colnames(x)
    an<-unique(c(cn,rn))
    myval<-x[!is.na(x)]
    mymat<-matrix(1,nrow=length(an),ncol=length(an),dimnames=list(an,an))
    for(ext in 1:length(cn))
    {
     for(int in 1:length(rn))
     {
     if(is.na(x[row.names(x)==rn[int],colnames(x)==cn[ext]])) next
     mymat[row.names(mymat)==rn[int],colnames(mymat)==cn[ext]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
     mymat[row.names(mymat)==cn[ext],colnames(mymat)==rn[int]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
     }
      
    }
    return(mymat)
}

if(to_do=="pca_analysis")
{
    run_pca(vcffile=vcffile,popfile=popfile)
} else if (to_do=="ibd_analysis")
{
    run_ibd(gdsfile=gdsfile,coefficient=coefficient,popfile=popfile,between=between,within=within,ibdmatfile=ibdmatfile,maf=maf,missingness=missingness,LD=LD)
}



