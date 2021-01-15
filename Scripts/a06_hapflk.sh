# Copyright:	    Gabriele Magris & Fabio Marroni 2020
# Aim:              Perform hapflk analysis 
# To add:		     
# Suggestions: 
# Fixes:  

source variable_list.txt

mkdir -p ${INPUT_DIR}/salmo_trutta_id1492/hapflk/output/
cd ${INPUT_DIR}/salmo_trutta_id1492/

# set parameters
K=20
nfit=20

# a. Run hapflk 
${HAPFLK_PATH}/hapflk/hapflk --file admixture/input/populations.snps.cov5.info50.sorted.no_low_cov -K ${K} --nfit=${nfit} --ncpu=4 --eigen --reynolds -p hapflk/output/populations.snps.cov5.info50.sorted.no_low_cov.K${K}.nfit${nfit} --kfrq &

# b. calculate p-value from hapflk values 
python ${HAPFLK_PATH}/utils/scaling_chi2_hapflk.py hapflk/output/populations.snps.cov5.info50.sorted.no_low_cov.K${K}.nfit${nfit}.hapflk ${K} 5

# c. draw the genome wide plot 
Rscript ${FUNC_DIR}/a06_hapflk.r \
-o ${INPUT_DIR}/salmo_trutta_id1492/hapflk2/plot/ \
-i ${INPUT_DIR}/salmo_trutta_id1492/hapflk2/output/populations.snps.cov5.info50.sorted.no_low_cov.K${K}.nfit${nfit}.hapflk_sc \
-c ${REFERENCE_DIR}/GCA_901001165.1_fSalTru1.1_genomic.fna.fai

