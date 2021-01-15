#!/bin/bash
#These are the main parameters you need to change to correctly set paths
FUNC_DIR=/folder/in/which/R/functions/are/located
INPUT_DIR=/base/directory/of/the/project/


#Read Zhp output and isolate results above threshold (Default threshold is 2.81 i.e. 0.0025 per tail, i.e. 0.005 in total, quite conservative, but not too much)

for EST in D f
do
Rscript $FUNC_DIR/Gensal/b13_ABBA_jack.r \
-I $INPUT_DIR/salmo_trutta_id1492/patterns_of_divergence/output/ \
-E $EST \
-O $INPUT_DIR/salmo_trutta_id1492/tables/ABBA_boot_${EST}.txt
done


#As performed for Edo, we measure "introgression" in the samples and add this information to Table_S1
#I did a ugly thing, I overwrote the file! Maybe you want to be cleaner than me!
Rscript $FUNC_DIR/Gensal/b01_extract_structure_4_edo.r \
-X $INPUT_DIR/salmo_trutta_id1492/tables/Table_S1.xlsx \
-A $INPUT_DIR/salmo_trutta_id1492/admixture/output/populations.snps.cov5.info50.sorted.no_low_cov.5.Qordered.population.origin.txt \
-O $INPUT_DIR/salmo_trutta_id1492/tables/Table_S1.xlsx

#Compute some basic statistics on SNP number and density
VCF=$INPUT_DIR/salmo_trutta_id1492/stacks/populations.snps.cov5.info50.no_low_cov.reheader.vcf
cut -f1 $VCF | sort | uniq -c > $INPUT_DIR/salmo_trutta_id1492/tables/SNP_stats.txt