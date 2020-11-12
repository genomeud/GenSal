# Run with --help or -h flag for help.
# Written 04/07/2020 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--infile"), type="character", default="",
  help="Formatted GO+KEGG file [default= %default]", metavar="character"), 
  make_option(c("-O", "--outfile"), type="character", default="", 
  help="Formatted GO+KEGG+KEGG pathways file [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$infile)) {
  stop("WARNING: No infile specified with '-G' flag.")
} else {  cat ("Infile ", opt$infile, "\n")
  infile <- opt$infile  
  }

if (is.null(opt$outfile)) {
  stop("WARNING: No outfile specified with '-G' flag.")
} else {  cat ("Outfile ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }

assign.kegg<-function(infile,outfile)
{
library(data.table)
#library(KEGGprofile)
library(KEGGREST)
#Get all organisms present in KEGG and their code
#fullorg<-keggList("organism")
indata<-fread(infile,data.table=F)
indata$pKEGG_map<-indata$pKEGG_class<-""
#We loop over the file with all the annotated transcript names (i.e. the names of transcript which returned a hit in blastx)
for(aaa in 1:nrow(indata))
{
	cat("Record", aaa, "out of",nrow(indata),"\n")	
	#Only analyze transcripts that have been associated to a KEGG term in the previous analysis
	#(note: KEGG term is in the form org:term)
	if(indata$KEGG[aaa]=="") next
	#Eventually split strings composed by more than one KEGG term (I think they are very rare, but just in case)
	inlist<-unlist(strsplit(indata$KEGG[aaa],";"))
	kmap_fin<-kclass_fin<-rep("",length(inlist))
	for(bbb in 1:length(inlist))
	{
	#See if we can associate a pathway to the KEGG term
	psname<-unique(keggLink("pathway", inlist[bbb]))
	if(length(psname)<1) next
	kmap_int<-kclass_int<-rep("",length(psname))
	for(ccc in 1:length(psname))
	{
		pfull<-keggGet(psname[ccc])
		#Read all the info on the full pathway (pfull)
		#And extract the info regarding the MAP and the CLASS, that can be used in the future.
		if(length(pfull)<1) next
		kclass<-kmap<-rep("",length(pfull))
		for(ddd in 1:length(pfull))
		{
		if(!is.null(pfull[[ddd]]$PATHWAY_MAP)) kmap[ddd]<-pfull[[ddd]]$PATHWAY_MAP
		if(!is.null(pfull[[ddd]]$CLASS)) kclass[ddd]<-pfull[[ddd]]$CLASS
		}
		kmap_int[ccc]<-paste(unique(kmap),collapse="|")
		kclass_int[ccc]<-paste(unique(kclass),collapse="|")
	}
	kmap_fin[bbb]<-paste(unique(kmap_int),collapse="|")
	kclass_fin[bbb]<-paste(unique(kclass_int),collapse="|")
	}
	#We use "|" as a separator because KEGG terms sometimes use ":" and ";" as separators.
	indata$pKEGG_map[aaa]<-paste(unique(kmap_fin),collapse="|")
	indata$pKEGG_class[aaa]<-paste(unique(kclass_fin),collapse="|")
}


write.table(indata,outfile,sep="\t",row.names=F,quote=F)

#To get pathway names of a KEGG term
#unique(keggLink("pathway","rno:84353"))

#To extract information of a pathway
#ff<-keggGet("path:rno04015")

#ciccio<-keggConv("uniprot","ath")
#Grep of the uniprot works!
#grep("Q8L6Y7",ciccio,value=T)

#This also works
#ciccio<-keggConv("ath","uniprot:Q8L6Y7")

#Superciccio is the best. We don't even need to know which one is the organism
#superciccio<-keggConv("genes", "uniprot:Q8L6Y7")

}
assign.kegg(infile=infile,outfile=outfile)
