# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
lfs quota -hg gwas_machine /path/to/project
# Disk quotas for grp gwas_machine (gid 15417):
#      Filesystem    used   quota   limit   grace   files   quota   limit   grace
# /path/to/project
#                      4k      0k      2T       -       1       0  101025       -

# singularity exec iibdgc_postprocess_10_singularity.sif

path="/path/to/project"

dataset=pQTL_decode

mkdir -p ${path}${dataset}/nominal/
mkdir -p ${path}${dataset}/permuted/
mkdir -p ${path}${dataset}/sample_size_per_condition/
mkdir -p ${path}${dataset}/log/
mkdir -p ${path}${dataset}/variant_info/

mkdir -p ${path}${dataset}/variant_info/raw/final_olink_ukb_bi/
mkdir -p ${path}${dataset}/variant_info/raw/final_somascan_smp/

cp /path/to/project ${path}${dataset}/variant_info/raw/final_olink_ukb_bi/
cp /path/to/project ${path}${dataset}/variant_info/raw/final_somascan_smp/

# recompress with bzgip and index with 
# bgzip
# tabix -f -s1 -b2 -e2 -S1

#######################################################################################################
# 0.- list complete list of studies in eQTL catalog:


#######################################################################################################
# 1.- reformat nominal summary stats file, and save in new location preserving the format Nikos uses:
# recreate path hierarchies by Nikos


# 1.1- Pararelize the creation of nominal file:

path_proteomics="/path/to/project"
path="/path/to/project"

MEM=2800

dataset=pQTL_decode
conditions=(final_olink_ukb_bi final_somascan_smp)

for condition in ${conditions[@]}
do 
cd ${path_proteomics}${condition}/ && files=($(ls | grep ".txt.bgz$" | sed 's/.txt.bgz//g'))  && echo  ${#files[@]}
for file in ${files[@]}
do
echo ${condition}_${file}
done
done


for condition in ${conditions[@]}
do 
cd ${path_proteomics}${condition}/ && files=($(ls | grep ".txt.bgz$" | sed 's/.txt.bgz//g'))  && echo  ${#files[@]}
for file in ${files[@]}
do
bsub -J"eQTL" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path}pQTL_decode/log/${dataset}_${condition}_${file}_stderr \
-o ${path}pQTL_decode/log/${dataset}_${condition}_${file}_stdout \
"Rscript ~/git/snakemake_colocalisation/scripts/pQTL_decode/process_pQTL_decode_nominal_summary_stats.R ${dataset} ${condition} ${file} > \
${path}pQTL_decode/log/process_pQTL_decode_nominal_summary_stats_${dataset}_${condition}_${file}.Rout"
done
done

# 2931
# 5284

rm ~/tmp.txt
for condition in ${conditions[@]}
do 
cd ${path_proteomics}${condition}/ && files=($(ls | grep ".txt.bgz$" | sed 's/.txt.bgz//g')) && for file in ${files[@]}
do
echo ${file} >> ~/tmp.txt && tail -50 ${path}pQTL_decode/log/${dataset}_${condition}_${file}_stdout | grep -E "Exited|Successfully" >> ~/tmp.txt
done
done


cat ~/tmp.txt | grep "Exited" | wc -l
cat ~/tmp.txt | grep "Successfully" | wc -l

# due to long name needs to be manually submitted
# GBR_UKB_OLINK_OID21013_FUT3_FUT5_3_galactosyl_N_acetylglucosaminide_4_alpha_L_fucosyltransferase_FUT3_4_galactosyl_N_acetylglucosaminide_3_alpha_L_fucosyltransferase_FUT5_adjAgeSexBatPC_InvNorm_22122022

# concatenate all files per condition:

cd /path/to/project

zcat Proteomics_SMP_PC0_13011_20_RPS10_RS10_10032022.txt.gz | head -1 > header.txt

MEM=800
bsub -J"eQTL" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path}pQTL_decode/log/final_olink_ukb_bi_concatenate_stderr \
-o ${path}pQTL_decode/log/final_olink_ukb_bi_concatenate_stdout \
"cat header.txt <(cat GBR_UKB_* | zcat | sed '/^Chrom/d') | gzip > final_olink_ukb_bi.txt.gz"


zcat final_olink_ukb_bi.txt.gz | head
zcat final_olink_ukb_bi.txt.gz | wc -l



MEM=800
bsub -J"eQTL" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path}pQTL_decode/log/final_somascan_smp_concatenate_stderr \
-o ${path}pQTL_decode/log/final_somascan_smp_concatenate_stdout \
"cat header.txt <(cat Proteomics_SMP_PC0_* | zcat | sed '/^Chrom/d') | gzip > final_somascan_smp.txt.gz"

zcat final_somascan_smp.txt.gz | head


dataset=pQTL_decode
conditions=(final_olink_ukb_bi final_somascan_smp)

MEM=125000

for condition in ${conditions[@]}
do
bsub -J"eQTL" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q week \
-e ${path}pQTL_decode/log/${dataset}_${condition}_process_stderr \
-o ${path}pQTL_decode/log/${dataset}_${condition}_process_stdout \
"Rscript ~/git/snakemake_colocalisation/scripts/pQTL_decode/process_pQTL_decode_nominal_and_permuted_summary_stats.R ${dataset} ${condition} > \
${path}pQTL_decode/log/process_pQTL_decode_nominal_and_permuted_summary_stats_${dataset}_${condition}.Rout"
done


##############################################
# 2.- generate sample sizes:


dataset=pQTL_decode
conditions=(final_olink_ukb_bi final_somascan_smp)

MEM=2000

for condition in ${conditions[@]}
do
bsub -J"eQTL" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path}pQTL_decode/log/${dataset}_${condition}_process_stderr \
-o ${path}pQTL_decode/log/${dataset}_${condition}_process_stdout \
"Rscript ~/git/snakemake_colocalisation/scripts/pQTL_decode/create_sample_size_file_pQTL_decode.R ${dataset} ${condition} > \
${path}pQTL_decode/log/create_sample_size_file_pQTL_decode_${dataset}_${condition}.Rout"
done

for condition in ${conditions[@]}
do
echo ${condition} && echo ${i} && tail -50  ${path}pQTL_decode/log/${dataset}_${condition}_process_stdout | grep -E "Successfully|Exited"
done


##############################################
# 3.- generate variant info files - completed in step 1


##############################################
# 5.- process the permuted files - already done in step 1


##############################################
# 6.-  final check that all files are in place:


for condition in ${conditions[@]}
do
echo ${condition} && \
ls -la ${path}${dataset}/nominal/${condition}/${condition}.tsv.gz && \
ls -la ${path}${dataset}/nominal/${condition}/${condition}.tsv.gz.tbi && \
ls -la ${path}${dataset}/permuted/${condition}.permuted.txt.gz && \
ls -la ${path}${dataset}/sample_size_per_condition/${condition} && \
ls -la ${path}${dataset}/variant_info/${condition}.variant.info.tsv.gz && \
ls -la ${path}${dataset}/variant_info/${condition}.variant.info.tsv.gz.tbi
done


##############################################
# 6.-  create a file listing the conditions:

for condition in ${conditions[@]}
do
echo ${condition} >>  ${path}${dataset}/list_conditions.txt
done

##############################################
# 7.-  submit coloc jobs 

export eqtl_study=${dataset}
export gwas_study=iibdgc

# TO BE EDITED"
# /path/to/project 
# eqtl_study variable NEEDS to be set up

mkdir -p /path/to/project

# MEM increased to 25K rather than 15K 
bsub -M 2000 -a "memlimit=True" -R "select[mem>2000] rusage[mem=2000] span[hosts=1]" \
-o /path/to/project \
-e /path/to/project \
-q oversubscribed -J "snakemake_master_COLOC" < /path/to/project

# evaluate if all jobs are completed correctly:
tail -50 /path/to/project | grep "Successfully"
head -10  /path/to/project
# job          count
# ---------  -------
# run_all          1
# run_coloc      120
# total          121

# re run
# job          count
# ---------  -------
# run_all          1
# run_coloc       60
# total           61

pheno=(IBD CD UC)

for ph in ${pheno[@]}
do
echo ${ph} && ls -la /path/to/project | wc -l 
done
# IBD
# 40
# CD
# 40
# UC
# 40

# evaluate integrity of files:
pheno=(IBD CD UC)


rm /path/to/project

for ph in ${pheno[@]}
do
echo ${ph} >> /path/to/project && \
for condition in ${conditions[@]}
do 
for chunck in {1..20}
do
echo ${condition} >> /path/to/project && \
head -1 /path/to/project | wc -w >> /path/to/project
tail -1 /path/to/project | wc -w >> /path/to/project
done
done
done

less /path/to/project

# combine all data:

export eqtl_study=${dataset}
export gwas_study=iibdgc

path="/path/to/project"
MEM=200

bsub -J"eQTL" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path}macromap/log/join_coloc_output.stderr \
-o ${path}macromap/log/join_coloc_output.stdout \
"bash /path/to/project"


unset eqtl_study
unset gwas_study

# empty log files folder:
rm -r /path/to/project
rm -r /path/to/project
