# Run with --help flag for help.
# Modified 08/11/2021 by Fabio Marroni

suppressPackageStartupMessages({
  library(optparse)
})


option_list = list(
  make_option(c("-I", "--infile"), type="character", default="",
              help="Input file", metavar="character"),
  make_option(c("-O", "--outfile"), type="character", default="", 
              help="output file name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

  if (is.null(opt$infile)) {
  stop("WARNING: No infile specified with '-I' flag.")
} else {  cat ("infile is ", opt$infile, "\n")
  infile <- opt$infile  
  }

  if (is.null(opt$outfile)) {
  stop("WARNING: No outfile specified with '-O' flag.")
} else {  cat ("outfile is ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }

plot_IT<-function(infile,outfile)
{
library("data.table")
library(mapIT)
library(ggplot2)
library(dplyr)
require(maps)
library(biogeo)
require(viridis)

theme_set(
  theme_void()
  )
library(openxlsx)
myloc<-read.xlsx(infile)
#Temp stuff to compute longitude in decimal degrees from GPS minutes and seconds
myloc$l<-unlist(lapply(strsplit(myloc$GPS.coordinates," "),"[",2))
myloc$ld<-as.numeric(unlist(lapply(strsplit(myloc$l,"°"),"[",1)))
myloc$lm<-as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(myloc$l,"°"),"[",2)),"'"),"[",1)))
myloc$ls<-as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(myloc$l,"'"),"[",2)),"\""),"[",1)))
myloc$lc<-substr(myloc$l,nchar(myloc$l),nchar(myloc$l))
#Longitude
myloc$long<-dms2dd(myloc$ld,myloc$lm,myloc$ls,myloc$lc)
#Temp stuff to compute latitude in decimal degrees from GPS minutes and seconds
myloc$l<-unlist(lapply(strsplit(myloc$GPS.coordinates," "),"[",1))
myloc$ld<-as.numeric(unlist(lapply(strsplit(myloc$l,"°"),"[",1)))
myloc$lm<-as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(myloc$l,"°"),"[",2)),"'"),"[",1)))
myloc$ls<-as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(myloc$l,"'"),"[",2)),"\""),"[",1)))
myloc$lc<-substr(myloc$l,nchar(myloc$l),nchar(myloc$l))
#Latitude
myloc$lat<-dms2dd(myloc$ld,myloc$lm,myloc$ls,myloc$lc)
mypoints<-data.frame(long=myloc$long,lat=myloc$lat,pop=myloc$Code)
mypoints<-mypoints[!duplicated(mypoints),]
mypoints$Lineage<-"orchid2"
mypoints$Lineage[mypoints$pop=="MA"]<-"gray68"
mypoints$Lineage[mypoints$pop=="AT"]<-"green1"
mypoints$Lineage[mypoints$pop=="MI"]<-"blue"
mypoints$Lineage[mypoints$pop=="MM"]<-"orangered"
mypoints$group<-NA


ciccio<-map_data("world")
my_map<-ciccio[ciccio$region=="Italy"|ciccio$subregion=="Corsica"|ciccio$region=="Austria",]
it_map <- map_data("world",region="Italy")
it_map<-my_map
#Add GPS coordinates of my samples to the map of Italy

mymap<- ggplot(it_map, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill="white", colour = "black")
fullmap<- mymap + geom_point(data=mypoints,aes(x=long, y=lat,color=Lineage), size=4,alpha=8/10) + 
  scale_color_manual(breaks=c("orchid2","gray68","green1","blue","orangered"),values=c("orchid2","gray68","green1","blue","orangered"),
  labels=c("CA","MA","AT","MI","MM")) +
  scale_shape_manual(breaks=c("orchid2","gray68","green1","blue","orangered"),values=c("orchid2","gray68","green1","blue","orangered"),
  labels=c("CA","MA","AT","MI","MM")) + 
  guides (size=F) #Removes size legend
ggsave(outfile,fullmap,device="png",type="cairo")
browser()
}

plot_IT(infile=infile,outfile=outfile)
