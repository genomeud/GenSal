# Copyright:	    Gabriele Magris & Fabio Marroni 2020
# Aim:              Perform phylogenetic analysis 
# To add:		     
# Suggestions: 
# Fixes:  

source variable_list.txt

mkdir -p ${INPUT_DIR}/salmo_trutta_id1492/snphylo/input/
cd ${INPUT_DIR}/salmo_trutta_id1492/

# Convert the chromosome names 
# create a conversion matrix with the chromosome name replaced by a number  
grep -v "^#" ${INPUT_DIR}/salmo_trutta_id1492/stacks/populations.snps.cov5.info50.no_low_cov.vcf  | \
cut -f 1 | \
uniq | \
# insert a number in the last column 
awk 'OFS="\t" {$(NF+1)=++i;}1' | \
# use this as new dictionary for the vcf translation 
awk  'FNR==NR{a[$1]=$2;next}{print a[$1],$0}'  - <(grep -v "^#" ${INPUT_DIR}/salmo_trutta_id1492/stacks/populations.snps.cov5.info50.no_low_cov.vcf) | \
sed 's/ /\t/1' | \
# remove the previous CHR ID 
cut -f 2 --complement | \
# append back the header 
cat <(grep "^#" ${INPUT_DIR}/salmo_trutta_id1492/stacks/populations.snps.cov5.info50.no_low_cov.vcf) - > snphylo/input/populations.snps.cov5.info50.no_low_cov.new_chr_names.vcf

# option -a the last autosome chr number 
cd snphylo/output/
# /projects/novabreed/share/gmagris/software/SNPhylo/snphylo.sh -v snphylo/input/populations.snps.new_chr_names.vcf -r -M 0.9 -m 0.01 -l 1 -P salmo_trutta_id1492.tree -b 100 -a 310

${SNPhylo_path}/snphylo.sh -v ../../snphylo/input/populations.snps.cov5.info50.no_low_cov.new_chr_names.vcf -r -M 0.8 -m 0.05 -l 1 -P salmo_trutta_id1492.cov5.info50.no_low_cov.tree -b 100 -a 310

# translate the names to article names 
prefix_output=salmo_trutta_id1492.cov5.info50.no_low_cov.tree
cp snphylo/output/${prefix_output}.ml.tree  snphylo/output/${prefix_output}.ml.bk.tree 
sed -i "s/:/|:/g" snphylo/output/${prefix_output}.ml.tree 
python /projects/novabreed/SNP/gatk/all_tmp_4/scripts/translate_id.py snphylo/output/translation_dictionary.txt snphylo/output/${prefix_output}.ml.tree


# perform bootstrap analysis (in R)
Rscript ${FUNC_DIR}/a03_SNPhylo.r \
--inpath ${INPUT_DIR}/salmo_trutta_id1492/snphylo/output/ \
--prefix salmo_trutta_id1492.cov5.info50.no_low_cov.tree
