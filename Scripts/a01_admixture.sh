# Copyright:	    Gabriele Magris & Fabio Marroni 2020
# Aim:              Perform admixture analysis + draw plot 
# To add:		     
# Suggestions: 
# Fixes:  

source variable_list.txt

mkdir -p ${INPUT_DIR}/salmo_trutta_id1492/admixture/input
cd ${INPUT_DIR}/salmo_trutta_id1492/admixture

# ----------------- #
# Create input file #
# ----------------- #

# need to transform the ped file obtained from stacks - first create bed file, because of an error in plink 
${PLINK_path}/plink --vcf ${INPUT_DIR}/salmo_trutta_id1492/stacks/populations.snps.cov5.info50.no_low_cov.vcf  --maf 0.000000001 --allow-extra-chr --out input/populations.snps.cov5.info50.no_low_cov --make-bed

# transform the bed file to ped /map 
${PLINK_path}/plink --bfile input/populations.snps.cov5.info50.no_low_cov --allow-extra-chr --recode 12 --out input/populations.snps.cov5.info50.no_low_cov

# sort the ped file in order to have the data separately for each population 
# assign to each individual the population name, sort the individuals and run admixture <- first column can be replaced with the family ID 
awk 'FNR==NR{a[$1]=$5;next}{ print $0, a[$1]}' <( sed 's/ /|/g' ${INPUT_DIR}/salmo_trutta_id1492/salmo_trutta_id1492_population_assignment.5groups.txt) input/populations.snps.cov5.info50.no_low_cov.ped  | \
# replace first column with last column and remove 
awk '{print $NF,$0}' | cut -f 2 --complement -d " " | \
# print range of columns 
awk '{out=$1;for(i=2;i<=(NF-1);i++){out=out" "$i}; print out}' | \
# sort the dataframe by population name 
sort -V -k1,1 > input/populations.snps.cov5.info50.no_low_cov.sorted.ped


# --------- #
# Admixture #
# --------- #

# 1. run Admixture 
mkdir -p ${INPUT_DIR}/salmo_trutta_id1492/admixture2/output/
cd ${INPUT_DIR}/salmo_trutta_id1492/admixture2/output/
for K in {1..10}
do
	${ADMIXTURE_path}/admixture -j24  ../input/populations.snps.cov5.info50.no_low_cov.sorted.ped ${K}
done


# 2. draw admixture plot 
Rscript ${FUNC_DIR}/a01_admixture.r \
-i ${INPUT_DIR}/salmo_trutta_id1492/admixture/input/ \
-o ${INPUT_DIR}/salmo_trutta_id1492/admixture/output/ \
--prefix populations.snps.cov5.info50.no_low_cov.sorted \
-k 5.Q \
-P ${INPUT_DIR}/salmo_trutta_id1492/salmo_trutta_id1492_population_assignment.5groups.txt \
--step draw_admixture


# 3. perform Evanno analysis + Cross-validation
mkdir -p evanno/

for a in {1..20}
do
    cd ${INPUT_DIR}/salmo_trutta_id1492/admixture/output/evanno
    mkdir -p ${a}
    cd ${a}
    for K in {1..16}
    do
        b=$((RANDOM*200000+10))
        ${ADMIXTURE_path}/admixture -s ${b} --cv=10 -j24 ../../../input/populations.snps.cov5.info50.sorted.no_low_cov.ped $K | tee log${K}.${a}.cov5.info50.no_low_cov.out
    done
done

# 3. collect results 
cd ${INPUT_DIR}/salmo_trutta_id1492/admixture/output/evanno
for run in {1..20}
do
    rm ${run}_cv-Loglikelihood.cov5.info50.no_low_cov.txt
    for k in {1..15}
    do
        likelihood=`grep "Loglikelihood" ${run}/log${k}.${run}.cov5.info50.no_low_cov.out | tail -n 1 | cut -f 2 -d " "`
        cv=`grep CV ${run}/log${k}.${run}.cov5.info50.no_low_cov.out | cut -f 4 -d " "`
        echo -e $k"\t"$cv"\t"$likelihood >> ${run}_cv-Loglikelihood.cov5.info50.no_low_cov.txt
    done
done

# 5. Draw CV plot - from the summary table
Rscript ${FUNC_DIR}/a01_admixture.r \
--infile ${INPUT_DIR}/salmo_trutta_id1492/admixture/output/populations.snps.sorted.no_low_cov.evanno.txt \
-o ${INPUT_DIR}/salmo_trutta_id1492/admixture3/output/ \
--prefix populations.snps.cov5.info50.no_low_cov.sorted \
--step draw_cross-validation


