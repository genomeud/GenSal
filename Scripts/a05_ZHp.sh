# Copyright:	    Gabriele Magris & Fabio Marroni 2020
# Aim:              Perform ZHp equation + draw genome-wide manhattan plot 
# To add:		     
# Suggestions: 
# Fixes:  

source variable_list.txt

mkdir -p ${INPUT_DIR}/salmo_trutta_id1492/selective_sweep/input
infile=${INPUT_DIR}/salmo_trutta_id1492/stacks/populations.snps.cov5.info50.no_low_cov.vcf

# extract for each sample the allele counts for reference/alternative
cd ${INPUT_DIR}/salmo_trutta_id1492/
for var in all_populations_cov5_info50_no_low_cov
do
    # # # cols=($(sed '1!d;s/\t/\n/g' <(head -n 15  ${infile}  | grep "^#CHR") | grep -nf <(grep ${var} salmo_trutta_id1492_population_assignment.5groups.txt |cut -f 1) | sed 's/:.*$//'))
    cols=($(sed '1!d;s/\t/\n/g' <(head -n 15  ${infile}  | grep "^#CHR") | grep -nf <(cut -f 1 salmo_trutta_id1492_population_assignment.5groups.txt | sed 1d ) | sed 's/:.*$//'))
    mkdir -p selective_sweep/input/${var}
    # loop now across each sample in the population that is stored in the array 
    for i in ${cols[@]}
    do
        sample_name=$(sed '1!d;s/\t/\n/g' <(head -n 15  ${infile}  | grep "^#CHR") | head -n ${i} | tail -n 1)
        # extract only the column of interest and split the field to get the allele counts -> REFERENCE
        cut -f 1,2,${i} ${infile} | \
        # remove the header (not necessary)
        grep -v "^#" | \
        # split the fields of interest 
        awk 'OFS="\t" {split($NF,fields,":");split(fields[3],allele,","); print allele[1]}' > selective_sweep/input/${var}/${sample_name}.reference.cov5.info50.no_low_cov.txt
        
        # extract only the column of interest and split the field to get the allele counts -> ALTERNATIVE 
        cut -f 1,2,${i} ${infile} | \
        # remove the header (not necessary)
        grep -v "^#" | \
        # split the fields of interest 
        awk 'OFS="\t" {split($NF,fields,":");split(fields[3],allele,","); print allele[2]}' > selective_sweep/input/${var}/${sample_name}.alternative.cov5.info50.no_low_cov.txt
        echo ${sample_name}
    done
echo ${var}
done


# I have created now a file for each sample with the counts of the reads supporting the reference or alternative allele.
# merge all populations and sum the reads counts per allele 
var=all_populations_cov5_info50_no_low_cov
out_var=all_populations
paste selective_sweep/input/${var}/*.reference.cov5.info50.no_low_cov.txt | awk '{sum=0; for (i=1; i<=NF; i++) { sum+= $i } print sum}' > selective_sweep/output/${out_var}.reference.reads.count.cov5.info50.no_low_cov.txt
paste selective_sweep/input/${var}/*.alternative.cov5.info50.no_low_cov.txt | awk '{sum=0; for (i=1; i<=NF; i++) { sum+= $i } print sum}' > selective_sweep/output/${out_var}.alternative.reads.count.cov5.info50.no_low_cov.txt

# now that the sum has been calculated -> reference + allele may be merged 
# need to add the coordinates in order to perform the windows analysis 
# Equation for Hp is 2*sum(Nmaj)*sum(Nmin) / (sum(Nmaj)+sum(Nmin))^2 
# Then the values is Z transformed ---> SCALED value (with average = 0, sd = 1) -> obtained: Hp-ave(Hp)/sd(Hp)
# in R it is obtained as :  scale(x, center = TRUE, scale = TRUE)

# for each population create a bed file with the coordinates and the reference / alternative counts 
paste <(grep -v "^#" ${infile} | awk 'OFS="\t" {print $1,$2-1,$2}') selective_sweep/output/${out_var}.reference.reads.count.cov5.info50.no_low_cov.txt selective_sweep/output/${out_var}.alternative.reads.count.cov5.info50.no_low_cov.txt | \
# append the header 
cat <(echo -e "chr\tstart\tend\treference\talternative") - > selective_sweep/output/${out_var}.merged.reads.count.cov5.info50.no_low_cov.txt
rm selective_sweep/output/${out_var}.reference.reads.count.cov5.info50.no_low_cov.txt selective_sweep/output/${out_var}.alternative.reads.count.cov5.info50.no_low_cov.txt
echo ${out_var}

# create the windows sizes and merge the data 
reference=${REFERENCE_DIR}/GCA_901001165.1_fSalTru1.1_genomic.fna
win_size=1000000
sliding=200000

# get the chromosomes size from the fai file and create the windows files 
bedtools makewindows -w ${win_size} -s ${sliding} -g <(cut -f 1,2 ${reference}.fai) -i srcwinnum | \
# remember that windows are sliding therefore bedtools merge cannot be used 
awk 'OFS="\t" {print $1,$2+1,$3,$4}' | \
# merge with the populations counts 
bedtools intersect  -wao -a - -b selective_sweep/output/${out_var}.merged.reads.count.cov5.info50.no_low_cov.txt | \
# sum based on same ID 
awk 'OFS="\t" {sum_ref[$1"|"$2"|"$3"|"$4] += $(NF-2); sum_alt[$1"|"$2"|"$3"|"$4] += $(NF-1)}; END{ for (id in sum_ref) { print id, sum_ref[id],sum_alt[id]} }' | \
sed 's/|/\t/g' | \
sort -k1,1 | \
# sort the df <- use the fai file to sort the chr <- because they are sorted by size 
join -a 1 <(cut -f 1 ${reference}.fai | awk 'OFS="\t" {$(NF+1)=++i;}1'  | grep "LR5844" | sort -k1,1 ) - | \
# replace spaces 
sed 's/ /\t/g' | \
# sort back in the correct order the df 
sort -V -k2,2 -k3n | \
# sort the columns 
awk 'OFS="\t" {print $1,$3,$4,$2,$5,$6,$7}' | \
# append the header
cat <(echo -e "chr\tstart\tend\tsort_order\twin_name\treference\talternative") - > selective_sweep/output/windows.${win_size}.sliding.${sliding}.${out_var}.merged.reads.count.cov5.info50.no_low_cov.txt
echo ${out_var}


# peform the ZHp equation 
Rscript ${FUNC_DIR}/Gensal/a05_ZHp.r \
-d ${INPUT_DIR}/selective_sweep/output/ \
-i ${INPUT_DIR}/selective_sweep/output/windows.${win_size}.sliding.${sliding}.${out_var}.merged.reads.count.cov5.info50.no_low_cov.txt \
-s ZHp_equation

# draw the genome wide plot 
Rscript ${FUNC_DIR}/Gensal/a05_ZHp.r \
-d ${INPUT_DIR}/selective_sweep/output/ \
-i ${INPUT_DIR}/selective_sweep/output/windows.${win_size}.sliding.${sliding}.${out_var}.merged.reads.count.cov5.info50.no_low_cov.ZHp.txt \
-s draw_plot
