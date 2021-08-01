#!/bin/bash
#To use these scripts you may need to play a bit with paths and folders
#
#These are the main parameters you need to change to correctly set paths
#Specify folders containing R functions and shell scripts downloaded from github 
FUNC_DIR=/folder/in/which/R/functions/are/located
SCRIPT_DIR=/folder/in/which/shell/scripts/are/located
INPUT_DIR=/base/directory/of/the/project/
BLASTDIR=/parent/folder/of/swissprot/blast/database
GODIR=/folder/in/which/GO/gaf/file/will/be/downloaded



##############################################################
#You will need to download Table_S1.xlsx of our paper to perform this analysis
#As performed for Edo, we measure "introgression" in the samples and add this information to Table_S1
#I did a ugly thing, I overwrote the file! Maybe you want to be cleaner than me!
Rscript $FUNC_DIR/Gensal/b01_extract_structure_4_edo.r \
-X $INPUT_DIR/salmo_trutta_id1492/tables/Table_S1.xlsx \
-A $INPUT_DIR/salmo_trutta_id1492/admixture/output/populations.snps.cov5.info50.sorted.no_low_cov.5.Qordered.population.origin.txt \
-O $INPUT_DIR/salmo_trutta_id1492/tables/Table_S1.xlsx
#################################################################


##############################################################
#You will need to download Table_S1.xlsx of our paper to perform this analysis
#Plot the map of Italy (plus Corsica and Austria) and puts the sampling points on it.
Rscript $FUNC_DIR/Gensal/b13_plot_map.r \
-I $INPUT_DIR/salmo_trutta_id1492/tables/Table_S1.xlsx \
-O $INPUT_DIR/salmo_trutta_id1492/tables/Fig_S1_Sampling.png
#################################################################




#Compute some basic statistics on SNP number and density
VCF=$INPUT_DIR/salmo_trutta_id1492/stacks/populations.snps.cov5.info50.no_low_cov.reheader.vcf
cut -f1 $VCF | sort | uniq -c > $INPUT_DIR/salmo_trutta_id1492/tables/SNP_stats.txt


###########################
#Main ZHp analysis
###########################
#Read Zhp output and isolate results above threshold (Default threshold is 2.81 i.e. 0.0025 per tail, i.e. 0.005 in total, quite conservative, but not too much)
Rscript $FUNC_DIR/Gensal/b02_top_ZHp.r \
-Z $INPUT_DIR/salmo_trutta_id1492/selective_sweep/output/windows.1000000.sliding.200000.all_populations.merged.reads.count.cov5.info50.no_low_cov.ZHp.txt \
-O $INPUT_DIR/salmo_trutta_id1492/selective_sweep/output/genes_windows.1000000.sliding.200000.all_populations.merged.reads.count.cov5.info50.no_low_cov.ZHp.txt \
-C $INPUT_DIR/genome/GCA_GCF_conversion.txt \
-G $INPUT_DIR/genome/GCF_901001165.1_fSalTru1.1_genomic.gff
#Compute significance of clustering of ZHp windows
Rscript $FUNC_DIR/Gensal/b15_sim_dist_ZHp.r \
-Z $INPUT_DIR/salmo_trutta_id1492/selective_sweep/output/windows.1000000.sliding.200000.all_populations.merged.reads.count.cov5.info50.no_low_cov.ZHp.txt \
-T 2.81 \
-S 1000 \
-L 27 \
-O $INPUT_DIR/salmo_trutta_id1492/tables/sim.ZHp.txt


#Download and gunzip the Uniprot GO:
cd ${GODIR}
wget -c ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz
gunzip goa_uniprot_all.gaf.gz


#Run blast on swissprot

# I decided to use the rna from genomics file (the largest, so I assume it has more features) to annotate.
# I checked that it has the same name convention of the gtf (obviously!!)
# 1. blastx against uniprot
cd $SCRIPT_DIR/logs
BLASTDB=$BLASTDIR/swissprot
MYQ=$INPUT_DIR/genome/GCF_901001165.1_fSalTru1.1_rna.fna
MYB=$INPUT_DIR/GO/GCF_901001165.1_fSalTru1.1_rna.out
MYGO=$GODIR/goa_uniprot_all.gaf
#module load aligners/blast/latest is installation specific. If you have blastx in your path, you can remove this part
#Also, qsub is needed to submit to a torque like queue system. You may need to remove this part
blast1=`echo "module load aligners/blast/latest; blastx -db $BLASTDB -query $MYQ -num_threads 12\
 -max_target_seqs 1 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' -evalue 1e-10 -out $MYB" | \
qsub -N b_rep -l vmem=20G,walltime=168:00:00,nodes=1:ppn=12`

##########################################################################
#Extract the Swissprot ID from blast output and write to file
##########################################################################
BSWISS=$INPUT_DIR/GO/blasted_swissprot_ID.txt
cut -f2 $MYB | cut -d"|" -f4| cut -d"." -f1 | sort | uniq > $BSWISS

##########################################################################
#Read the gene_id in $BSWISS and return them if they have a correspondence in the second column of $MYGO
##########################################################################

MYMAPPING=$INPUT_DIR/GO/transcripts_GO.txt
#This might be very slow, and we might rather run the command below on the split files.
awk 'FNR==NR{k[$1]=1;next;} k[$2]' $BSWISS $MYGO > $MYMAPPING 


#Reformat GO file, so that each gene has only one row and all GO terms are separated by ";"
#We recycle a function written for the work on limodorum
FORMGO=$INPUT_DIR/GO/formatted_transcripts_GO.txt
Rscript $FUNC_DIR/Gensal/b03_reformat_GO.r \
-I $MYMAPPING \
-O $FORMGO

#Reformat GO file, so that each gene has only one row and all GO terms are separated by ";"
#We recycle a function written for the work on limodorum
KEGGO=$INPUT_DIR/GO/formatted_transcripts_GO_KEGG.txt
Rscript $FUNC_DIR/Gensal/b04_assign_kegg.r \
-I $FORMGO \
-O $KEGGO

#Reformat KEGG file, so that each gene has only one row and all GO terms are separated by ";"
#We recycle a function written for the work on limodorum
PKEGG=$INPUT_DIR/GO/formatted_transcripts_GO_KEGG_path.txt
Rscript $FUNC_DIR/Gensal/b05_assign_kegg_pathway.r \
-I $KEGGO \
-O $PKEGG

#Add original transcript name (useful to associate my genes to KEGG and GO!!!)
WITHT=$INPUT_DIR/GO/named_transcripts_GO_KEGG_path.txt
Rscript $FUNC_DIR/Gensal/b06_assign_transcript_name.r \
-I $PKEGG \
-B $MYB \
-O $WITHT

#Perform KEGG enrichment on ZHp results
ENRICH=$INPUT_DIR/GO/ZHp_KEGG_Enrich.txt
Rscript $FUNC_DIR/Gensal/b07_KEGG_enrichment.r \
-I $INPUT_DIR/salmo_trutta_id1492/selective_sweep/output/genes_windows.1000000.sliding.200000.all_populations.no_low_cov.merged.reads.count.ZHp.txt \
-K $WITHT \
-O $ENRICH




###########################
#ZHp analysis on subgroups
###########################
#Read Zhp output and isolate results above threshold (Default threshold is 2.81 i.e. 0.0025 per tail, i.e. 0.005 in total, quite conservative, but not too much)
for GROUP in Q95 Farm River carpione fario_atlantica fario_med_island fario_med_peninsula marmorata
do
Rscript $FUNC_DIR/Gensal/b02_top_ZHp.r \
-Z $INPUT_DIR/salmo_trutta_id1492/selective_sweep/output/windows.1000000.sliding.200000.${GROUP}.merged.reads.count.cov5.info50.no_low_cov.ZHp.txt \
-O $INPUT_DIR/salmo_trutta_id1492/selective_sweep/output/genes_windows.1000000.sliding.200000.${GROUP}.merged.reads.count.cov5.info50.no_low_cov.ZHp.txt \
-C $INPUT_DIR/genome/GCA_GCF_conversion.txt \
-G $INPUT_DIR/genome/GCF_901001165.1_fSalTru1.1_genomic.gff
#Compute significance of clustering of ZHp windows
Rscript $FUNC_DIR/Gensal/b15_sim_dist_ZHp.r \
-Z $INPUT_DIR/salmo_trutta_id1492/selective_sweep/output/windows.1000000.sliding.200000.${GROUP}.merged.reads.count.cov5.info50.no_low_cov.ZHp.txt \
-T 2.81 \
-S 1000 \
-L 27 \
-O $INPUT_DIR/salmo_trutta_id1492/tables/sim.${GROUP}.ZHp.txt
done

for GROUP in Q95 Farm River carpione fario_atlantica fario_med_island fario_med_peninsula marmorata
do
#Perform KEGG enrichment on ZHp results
ENRICH=$INPUT_DIR/GO/ZHp_KEGG_Enrich_${GROUP}.txt
Rscript $FUNC_DIR/Gensal/b07_KEGG_enrichment.r \
-I $INPUT_DIR/salmo_trutta_id1492/selective_sweep/output/genes_windows.1000000.sliding.200000.${GROUP}.merged.reads.count.cov5.info50.no_low_cov.ZHp.txt \
-K $WITHT \
-O $ENRICH
done




#Check if we find some signal in mucins genes, on which the work on trout transcriptome focused(https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5315281/)
module load lang/perl/5.10.1
MUCINS=$INPUT_DIR/salmo_trutta_id1492/mucins/mucins.fasta
perl $SCRIPT_DIR/Gensal/use_eutils_mucins_onlyfasta.pl > $MUCINS
#At the end I didn't find anything interesting in the mucins and I didn't talk about this in the paper.


#Creat blastdb for salmo trutta
module load aligners/blast/latest
makeblastdb -in $INPUT_DIR/genome/GCF_901001165.1_fSalTru1.1_genomic.fna -dbtype nucl -out $INPUT_DIR/genome/GCF_901001165.1_fSalTru1.1_genomic
#Blast Mucins against salmo trutta genome
blastn -query $MUCINS -db $INPUT_DIR/genome/GCF_901001165.1_fSalTru1.1_genomic -evalue=10e-10 -max_target_seqs 5 -outfmt '6' -out ${MUCINS/.fasta/.out}


##############################
#
# HapFLK
#
##############################

#Extract significant results of HapFLK. We consider as signficant only windows consisting of at least 2 consecutive SNPs with pvalue lower than a threshold (0.01 by default).
#Windows are interrupted when a non-signficant SNP is met or when 2 consecutive significant SNPs are more than 1Mb apart.  
Rscript $FUNC_DIR/Gensal/b08_top_hapFLK.r \
-I $INPUT_DIR/salmo_trutta_id1492/hapflk/output/populations.snps.cov5.info50.sorted.no_low_cov.K20.nfit20.hapflk_sc \
-D 1000000 \
-P 0.01 \
-S 2 \
-C $INPUT_DIR/genome/GCA_GCF_conversion.txt \
-G $INPUT_DIR/genome/GCF_901001165.1_fSalTru1.1_genomic.gff \
-O $INPUT_DIR/salmo_trutta_id1492/tables/top_populations.snps.cov5.info50.sorted.no_low_cov.K20.nfit20.hapflk_sc


#Perform KEGG enrichment on HapFLK results
WITHT=$INPUT_DIR/GO/named_transcripts_GO_KEGG_path.txt
ENRICH=$INPUT_DIR/GO/HapFLK_KEGG_Enrich.txt
Rscript $FUNC_DIR/Gensal/b07_KEGG_enrichment.r \
-I $INPUT_DIR/salmo_trutta_id1492/tables/top_populations.snps.cov5.info50.sorted.no_low_cov.K20.nfit20.hapflk_sc \
-K $WITHT \
-O $ENRICH

#Check my distribution of clusters VS random distribution. First line is real data, all the other are random
#Done for hapFLK, need to do for ZHp. Also need to summarize the results
Rscript $FUNC_DIR/Gensal/b14_sim_dist_hapFLK.r \
-I $INPUT_DIR/salmo_trutta_id1492/hapflk/output/populations.snps.cov5.info50.sorted.no_low_cov.K20.nfit20.hapflk_sc \
-D 1000000 \
-P 0.01 \
-S 2 \
-L 50 \
-N 250 \
-O $INPUT_DIR/salmo_trutta_id1492/tables/sim.hapflk.txt

#####################################
#
#HapFLK on subgroups
#Farm vs river and Q95
#
#####################################


#Extract significant results of HapFLK. We consider as signficant only windows consisting of at least 2 consecutive SNPs with pvalue lower than a threshold (0.01 by default).
#Windows are interrupted when a non-signficant SNP is met or when 2 consecutive significant SNPs are more than 1Mb apart.  
for MYSET in farm_vs_river Q95
do
Rscript $FUNC_DIR/Gensal/b08_top_hapFLK.r \
-I $INPUT_DIR/salmo_trutta_id1492/hapflk/output/populations.snps.cov5.info50.sorted.no_low_cov.${MYSET}.subset.K20.nfit20.hapflk_sc \
-D 1000000 \
-P 0.01 \
-S 2 \
-C $INPUT_DIR/genome/GCA_GCF_conversion.txt \
-G $INPUT_DIR/genome/GCF_901001165.1_fSalTru1.1_genomic.gff \
-O $INPUT_DIR/salmo_trutta_id1492/tables/top_populations.snps.cov5.info50.sorted.no_low_cov.${MYSET}.subset.K20.nfit20.hapflk_sc
#Perform KEGG enrichment on HapFLK results
WITHT=$INPUT_DIR/GO/named_transcripts_GO_KEGG_path.txt
ENRICH=$INPUT_DIR/GO/HapFLK_${MYSET}_KEGG_Enrich.txt
Rscript $FUNC_DIR/Gensal/b07_KEGG_enrichment.r \
-I $INPUT_DIR/salmo_trutta_id1492/tables/top_populations.snps.cov5.info50.sorted.no_low_cov.${MYSET}.subset.K20.nfit20.hapflk_sc \
-K $WITHT \
-O $ENRICH
#Check my distribution of clusters VS random distribution. First line is real data, all the other are random
#Done for hapFLK, need to do for ZHp. Also need to summarize the results
Rscript $FUNC_DIR/Gensal/b14_sim_dist_hapFLK.r \
-I $INPUT_DIR/salmo_trutta_id1492/hapflk/output/populations.snps.cov5.info50.sorted.no_low_cov.${MYSET}.subset.K20.nfit20.hapflk_sc \
-D 1000000 \
-P 0.01 \
-S 2 \
-L 50 \
-N 250 \
-O $INPUT_DIR/salmo_trutta_id1492/tables/sim.${MYSET}.hapflk.txt
done







##############################
#
# Windows of IBD between populations
#
##############################

#Select the most interesting windows, i.e. those where IBD sharing between populations is greater
GFF=$INPUT_DIR/genome/GCF_901001165.1_fSalTru1.1_genomic.gff
IBDGENES=$INPUT_DIR/salmo_trutta_id1492/tables/genes_top_all_chr_kinship_median.coord.txt 
Rscript $FUNC_DIR/Gensal/b09_top_win_ibd.r \
-I $INPUT_DIR/salmo_trutta_id1492/SNPRelate/windows_ibd/kinship_median/all_chr_kinship_median.coord.txt \
-K 0.05 \
-C $INPUT_DIR/genome/GCA_GCF_conversion.txt \
-G $GFF \
-O $IBDGENES

#Use resampling to assess significance of the windows.
Rscript $FUNC_DIR/Gensal/b16_sim_dist_ibd.r \
-I $INPUT_DIR/salmo_trutta_id1492/SNPRelate/windows_ibd/kinship_median/all_chr_kinship_median.coord.txt \
-K 0.05 \
-L 27 \
-S 1000 \
-O $INPUT_DIR/salmo_trutta_id1492/tables/sim.ibd.txt

#Perform KEGG enrichment
WITHT=$INPUT_DIR/GO/named_transcripts_GO_KEGG_path.txt
ENRICH=$INPUT_DIR/GO/win_ibd_KEGG_Enrich.txt
Rscript $FUNC_DIR/Gensal/b07_KEGG_enrichment.r \
-I $IBDGENES \
-K $WITHT \
-O $ENRICH


#Compare farm vs river

#Select the most interesting windows, i.e. those where IBD sharing between populations is greater
GFF=$INPUT_DIR/genome/GCF_901001165.1_fSalTru1.1_genomic.gff
IBDGENES=$INPUT_DIR/salmo_trutta_id1492/tables/genes_top_all_chr_kinship_median.farm_vs_river.coord.txt 
Rscript $FUNC_DIR/Gensal/b09_top_win_ibd.r \
-I $INPUT_DIR/salmo_trutta_id1492/SNPRelate/windows_ibd/kinship_median/all_chr_kinship_median.farm_vs_river.coord.txt \
-K 0.05 \
-C $INPUT_DIR/genome/GCA_GCF_conversion.txt \
-G $GFF \
-O $IBDGENES
#Use resampling to assess significance of the windows.
Rscript $FUNC_DIR/Gensal/b16_sim_dist_ibd.r \
-I $INPUT_DIR/salmo_trutta_id1492/SNPRelate/windows_ibd/kinship_median/all_chr_kinship_median.farm_vs_river.coord.txt \
-K 0.05 \
-L 27 \
-S 1000 \
-O $INPUT_DIR/salmo_trutta_id1492/tables/sim.ibd.farm_vs_river.txt
#Perform KEGG enrichment
WITHT=$INPUT_DIR/GO/named_transcripts_GO_KEGG_path.txt
ENRICH=$INPUT_DIR/GO/win_ibd_farm_vs_river_KEGG_Enrich.txt
Rscript $FUNC_DIR/Gensal/b07_KEGG_enrichment.r \
-I $IBDGENES \
-K $WITHT \
-O $ENRICH

#Excluding admixed individuals

#Select the most interesting windows, i.e. those where IBD sharing between populations is greater
GFF=$INPUT_DIR/genome/GCF_901001165.1_fSalTru1.1_genomic.gff
IBDGENES=$INPUT_DIR/salmo_trutta_id1492/tables/genes_top_all_chr_kinship_median.Q95.coord.txt 
Rscript $FUNC_DIR/Gensal/b09_top_win_ibd.r \
-I $INPUT_DIR/salmo_trutta_id1492/SNPRelate/windows_ibd/kinship_median/all_chr_kinship_median.Q95.coord.txt \
-K 0.05 \
-C $INPUT_DIR/genome/GCA_GCF_conversion.txt \
-G $GFF \
-O $IBDGENES
#Use resampling to assess significance of the windows.
Rscript $FUNC_DIR/Gensal/b16_sim_dist_ibd.r \
-I $INPUT_DIR/salmo_trutta_id1492/SNPRelate/windows_ibd/kinship_median/all_chr_kinship_median.Q95.coord.txt \
-K 0.05 \
-L 27 \
-S 1000 \
-O $INPUT_DIR/salmo_trutta_id1492/tables/sim.ibd.Q95.txt
#Perform KEGG enrichment
WITHT=$INPUT_DIR/GO/named_transcripts_GO_KEGG_path.txt
ENRICH=$INPUT_DIR/GO/win_ibd_Q95_KEGG_Enrich.txt
Rscript $FUNC_DIR/Gensal/b07_KEGG_enrichment.r \
-I $IBDGENES \
-K $WITHT \
-O $ENRICH



#############################################
#
# Add chromosome numbers (I did it manually, but after having to redo that several times I finally decided to write a function)
#
#############################################

NAMETONUMBERS=$INPUT_DIR/genome/GCA_901001165.1_sequence_report.txt
ENRICH=$INPUT_DIR/GO/win_ibd_KEGG_Enrich.txt
Rscript $FUNC_DIR/Gensal/b10_add_Chr_Numb.r \
-N $NAMETONUMBERS \
-O $ENRICH

#Just for the sake of comparison I also computed Fis at the whole population level. Be careful Fis average will not be written to file but to stdout.
Rscript $FUNC_DIR/Gensal/b11_compute_Fis_He_Ho.r \
-H $INPUT_DIR/salmo_trutta_id1492/heterozygosity/populations.snps.cov5.info50.no_low_cov.hwe

#Read Zhp output and isolate results above threshold (Default threshold is 2.81 i.e. 0.0025 per tail, i.e. 0.005 in total, quite conservative, but not too much)
for EST in D f
do
Rscript $FUNC_DIR/Gensal/b12_ABBA_jack.r \
-I $INPUT_DIR/salmo_trutta_id1492/patterns_of_divergence/output/ \
-E $EST \
-O $INPUT_DIR/salmo_trutta_id1492/tables/ABBA_boot_${EST}.txt
done
