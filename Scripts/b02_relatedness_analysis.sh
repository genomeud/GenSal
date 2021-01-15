#!/bin/bash
#This are scripts for downstream analysis of the relatedness results.

#These are the main parameters you need to change to correctly set paths
FUNC_DIR=/folder/in/which/R/functions/are/located
INPUT_DIR=/base/directory/of/the/project/

#Use SNPRelate to compute PCA 
Rscript $FUNC_DIR/Gensal/b11_SNPrelate.r \
-V $INPUT_DIR/salmo_trutta_id1492/stacks/populations.snps.cov5.info50.no_low_cov.vcf \
-P $INPUT_DIR/salmo_trutta_id1492/salmo_trutta_id1492_population_assignment.5groups.txt \
# missing pca plot (jpeg)


#Use SNPRelate to compute IBD and use it as relatedness measure (I think I like it more than related coefficients)
#I tried several combinations to check for robustness and see if high relatedness was due to some artifact.
#At the end we were reassured that IBD was robust to changing analysis parameters, and we didn't show all these tests in the paper.
MAF=0.05
MISSING=0.05
LD=0.2

Rscript $FUNC_DIR/Gensal/b11_SNPrelate.r \
-G $INPUT_DIR/salmo_trutta_id1492/SNPRelate/salmo_trutta_id1492.cov5.info50.gds \
-P $INPUT_DIR/salmo_trutta_id1492/salmo_trutta_id1492_population_assignment.5groups.txt \
-B $INPUT_DIR/salmo_trutta_id1492/SNPRelate/IBD_between_maf_${MAF}_missing_${MISSING}_LD_${LD}.pdf \
-W $INPUT_DIR/salmo_trutta_id1492/SNPRelate/IBD_within_maf_${MAF}_missing_${MISSING}_LD_${LD}.pdf \
-H $INPUT_DIR/salmo_trutta_id1492/selective_sweep/input/all_populations_cov5_info50_no_low_cov/ \
-I $INPUT_DIR/salmo_trutta_id1492/SNPRelate/salmo_trutta_id1492.cov5.info50.ibd \
-M $MAF \
-S $MISSING \
-L $LD



#Just for the sake of comparison I also computed Fis at the whole population level. Be careful Fis average will not be written to file but to stdout.
Rscript $FUNC_DIR/Gensal/b12_compute_Fis_He_Ho.r \
-H $INPUT_DIR/salmo_trutta_id1492/heterozygosity/populations.snps.cov5.info50.no_low_cov.hwe

