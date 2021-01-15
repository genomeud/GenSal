# Copyright:	    Gabriele Magris & Fabio Marroni 2020
# Aim:              Perform treemix analysis 
# To add:		     
# Suggestions: 
# Fixes:  


source variable_list.txt

mkdir -p ${INPUT_DIR}/salmo_trutta_id1492/treemix2/input/
cd ${INPUT_DIR}/salmo_trutta_id1492/treemix2


# need to transform the ped file obtained from stacks - first create bed file, because of an error in plink 
${PLINK_path}/plink --vcf ${INPUT_DIR}/salmo_trutta_id1492/stacks/populations.snps.cov5.info50.no_low_cov.vcf  --maf 0.000000001 --allow-extra-chr --out input/populations.snps.cov5.info50.no_low_cov --make-bed

# transform the bed file to ped /map 
${PLINK_path}/plink --bfile input/populations.snps.cov5.info50.no_low_cov --allow-extra-chr --recode tab --out input/populations.snps.cov5.info50.no_low_cov



# transpose the ped file to get in each row a snp // each column a sample 
transpose <input/populations.snps.cov5.info50.no_low_cov.ped| \
# remove unwanted rows 
sed 1d | sed 2d | sed 2d | sed 2d | sed 2d | \
paste <(cat <(echo -e "NA\tNA\tNA\tNA") input/populations.snps.cov5.info50.no_low_cov.map) -  > input/populations.snps.cov5.info50.no_low_cov.transposed.txt

# create a file with the ID of the SNPs 
cut -f 2 input/populations.snps.cov5.info50.no_low_cov.transposed.txt  | sed 1d > input/populations.snps.cov5.info50.no_low_cov.all_chr.geno.selected.txt

group_assignement=${INPUT_DIR}/salmo_trutta_id1492/related/input/population_assignment.txt
mkdir -p all_chr/
for in_group in carpione fario_atlantica fario_med_island fario_med_peninsula marmorata
do
    # in_group=chip_BAL.txt
    outprefix=$(basename $in_group .txt)
    col=($(sed '1!d;s/\t/\n/g' <(head -n 1 input/populations.snps.cov5.info50.no_low_cov.transposed.txt ) | grep -nf <(grep -w ${in_group} ${group_assignement} | cut -f 1) | sed 's/:.*$//' | sed ':a;N;$!ba;s/\n/,/g'))

    # extract from the transposed file the column of interest 
    cut -f ${col}  input/populations.snps.cov5.info50.no_low_cov.transposed.txt  | \
    # calculate the allele frequency 
    awk '{print gsub(/A/,"")"\t"gsub(/C/,"")"\t"gsub(/G/,"")"\t"gsub(/T/,"")}' | \
    # remove header row
    sed 1d | \
    # append the coordinates from the file 
    paste <(cut -f 2  input/populations.snps.cov5.info50.no_low_cov.transposed.txt  | sed 1d ) - > all_chr/${outprefix}_allele_count.cov5.info50.txt

    # assign the correct orientation of the SNP <- use annovar for this purpose 
    # (because all coordinates need to be assigned in our dataset)
    paste <( cat  all_chr/${outprefix}_allele_count.cov5.info50.txt ) <(cut -f 4,5 ${INPUT_DIR}/salmo_trutta_id1492/stacks/populations.snps.cov5.info50.no_low_cov.vcf | grep -v "^#" | sed 1d ) | sed 's/ /\t/g' | awk 'BEGIN {OFS="\t"} {if(($6=="A" || $7=="A") && ($6=="C" || $7=="C")) print $0,$2,$3; if(($6=="A" || $7=="A") && ($6=="G" || $7=="G")) print $0,$2,$4; if(($6=="A" || $7=="A") && ($6=="T" || $7=="T")) print $0,$2,$5; if(($6=="C" || $7=="C") && ($6!="A" || $7!="A") && ($6=="G" || $7=="G")) print $0,$3,$4 ; if(($6=="C" || $7=="C") && ($6!="A" || $7!="A") && ($6=="T" || $7=="T")) print $0,$3,$5 ; if(($6=="G" || $7=="G") && ($6!="A" || $7!="A") && ($6!="C" || $7!="C") && ($6=="T" || $7=="T")) print $0,$4,$5  }' | awk 'BEGIN {OFS="\t"}{print $1,$8","$9}' > all_chr/${outprefix}_allele_count_ordered.cov5.info50.txt
    echo ${outprefix}
done

# merge now all necessary groups 
outprefix=K5.snp_chip.groups_allele_count_ordered.cov5.info50
tmp_path=all_chr
awk '{a[FNR] = a[FNR] (NR==FNR?"":"\t") $2} END { for (i=1;i<=FNR;i++) print a[i] }'  \
${tmp_path}/carpione_allele_count_ordered.cov5.info50.txt \
${tmp_path}/fario_atlantica_allele_count_ordered.cov5.info50.txt \
${tmp_path}/fario_med_island_allele_count_ordered.cov5.info50.txt \
${tmp_path}/fario_med_peninsula_allele_count_ordered.cov5.info50.txt \
${tmp_path}/marmorata_allele_count_ordered.cov5.info50.txt | \
paste input/populations.snps.cov5.info50.no_low_cov.all_chr.geno.selected.txt - | cut -f2- | grep -vP "\t0,0\t"  | grep -vP "^0,0\t" | grep -vP "\t0,0$" | sed 's/\t/ /g' >  ${tmp_path}/${outprefix}.txt

# create a file with groups order 
mkdir -p groups/
echo -e "Carpione\nAtlantica\nMediterranea|Island\nMediterranea|Mainland\nMarmoratus" > groups/poporder_K5_ddrad.txt

# add header
mkdir -p output
cat <(echo 'Carpione Atlantica Mediterranea|Island Mediterranea|Mainland Marmoratus' ) ${tmp_path}/${outprefix}.txt > output/${outprefix}.txt

gzip output/${outprefix}.txt


## RUN TREEMIX
cd output/

in_file=${outprefix}
out_file=K5_ddrad.cov5.info50

k=10
migration=3
${TREEMIX_path}/treemix -se -i ${in_file}.txt.gz -k ${k} -o ${out_file}_${k}_${migration}_se -m ${migration} --bootstrap &


# draw plot 
Rscript ${FUNC_DIR}/a04_treemix.r \
-d ${INPUT_DIR}/salmo_trutta_id1492/treemix2/output \
-p ${out_file}_${k}_${migration}_se \
-f ${FUNC_DIR}/



# run F3 F4 
in_file=${outprefix}
out_file=K5_ddrad.cov5.info50
k=10
${TREEMIX_path}/threepop -i ${in_file}.txt.gz -k ${k} > ${outprefix}.threepop.txt 
${TREEMIX_path}/fourpop -i ${in_file}.txt.gz -k ${k} > ${outprefix}.fourpop.txt 
