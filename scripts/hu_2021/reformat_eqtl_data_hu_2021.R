# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# singularity exec iibdgc_postprocess_10_singularity.sif

# download data and create new project directory:

path="/path/to/project"

mkdir -p ${path}QTL_managed_access/hu_2021/raw_data/
mkdir -p ${path}QTL_managed_access/hu_2021/other/
mkdir -p ${path}QTL_managed_access/hu_2021/log/
mkdir -p ${path}QTL_managed_access/hu_2021/nominal/
mkdir -p ${path}QTL_managed_access/hu_2021/permuted/
mkdir -p ${path}QTL_managed_access/hu_2021/sample_size_per_condition/
mkdir -p ${path}QTL_managed_access/hu_2021/variant_info/




#######################################################################################################
# 0.- dataset in b37, liftover to b38 first:

MEM=20000

bsub -J"eQTL_lift" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path}QTL_managed_access/hu_2021/log/hu_2021_liftover_process_stderr \
-o ${path}QTL_managed_access/hu_2021/log/hu_2021_liftover_process_stdout \
"Rscript ~/git/snakemake_colocalisation/scripts/hu_2021/liftover_hu_2021_b37_to_b38.R > \
${path}QTL_managed_access/hu_2021/log/liftover_hu_2021_b37_to_b38.Rout"

# Completed


##############################################
# 1.- Create sample size per trait files, variant info files, permuted and nominal files:

MEM=55000

bsub -J"eQTL_perm" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path}QTL_managed_access/hu_2021/log/hu_2021_nominal_permuted_process_stderr \
-o ${path}QTL_managed_access/hu_2021/log/hu_2021_nominal_permuted_process_stdout \
"Rscript ~/git/snakemake_colocalisation/scripts/hu_2021/process_hu_2021_nominal_and_permuted_summary_stats.R > \
${path}QTL_managed_access/hu_2021/log/process_hu_2021_nominal_and_permuted_summary_stats.Rout"


##############################################
# 2.-  final check that all files are in place:

cohort="intestinal_mucosa_colon_ileum"
  
echo ${cohort} && \
ls -la ${path}QTL_managed_access/hu_2021/nominal/${cohort}/${cohort}.tsv.gz && \
ls -la ${path}QTL_managed_access/hu_2021/nominal/${cohort}/${cohort}.tsv.gz.tbi && \
ls -la ${path}QTL_managed_access/hu_2021/permuted/${cohort}.permuted.txt.gz && \
ls -la ${path}QTL_managed_access/hu_2021/sample_size_per_condition/${cohort} && \
ls -la ${path}QTL_managed_access/hu_2021/variant_info/${cohort}.variant.info.tsv.gz && \
ls -la ${path}QTL_managed_access/hu_2021/variant_info/${cohort}.variant.info.tsv.gz.tbi

# intestinal_mucosa_colon_ileum
# -rw-rw---- 1 lf9 gwas_machine 697814645 Jun 23 09:29 /path/to/project
# -rw-rw---- 1 lf9 gwas_machine 1728314 Jun 23 09:30 /path/to/project
# -rw-rw---- 1 lf9 gwas_machine 1316385 Jun 23 09:34 /path/to/project
# -rw-rw---- 1 lf9 gwas_machine 34 Jun 23 09:35 /path/to/project
# -rw-rw---- 1 lf9 gwas_machine 34657487 Jun 23 09:35 /path/to/project
# -rw-rw---- 1 lf9 gwas_machine 1347901 Jun 23 09:35 /path/to/project



##############################################
# 6.-  create a file listing the conditions:

dataset="QTL_managed_access/hu_2021"

echo "intestinal_mucosa_colon_ileum" >>  ${path}${dataset}/list_conditions.txt


##############################################
# 7.-  submit coloc jobs 

export eqtl_study=hu_2021
export gwas_study=iibdgc

# TO BE EDITED"
# /path/to/project 
# eqtl_study variable NEEDS to be set up

# /path/to/project
# uncommnet line 21 for managed access projects

mkdir -p /path/to/project

bsub -M 2000 -a "memlimit=True" -R "select[mem>2000] rusage[mem=2000] span[hosts=1]" \
-o /path/to/project \
-e /path/to/project \
-q oversubscribed -J "snakemake_master_COLOC" < /path/to/project


# evaluate if all jobs are completed correctly:
tail -50 /path/to/project | grep "Successfully"
head -15  /path/to/project
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
# 20
# CD
# 20
# UC
# 20

# evaluate integrity of files:
pheno=(IBD CD UC)

rm /path/to/project

for ph in ${pheno[@]}
do
echo ${ph} >> /path/to/project && \
for chunck in {1..20}
do
head -1 /path/to/project | wc -w >> /path/to/project
tail -1 /path/to/project | wc -w >> /path/to/project
done
done

less /path/to/project

# combine all data:

export eqtl_study=hu_2021
export gwas_study=iibdgc

path="/path/to/project"
MEM=200

bsub -J"eQTL" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path}QTL_managed_access/${eqtl_study}/log/join_coloc_output.stderr \
-o ${path}QTL_managed_access/${eqtl_study}/log/join_coloc_output.stdout \
"bash /path/to/project"


rm /path/to/project
unset eqtl_study
unset gwas_study

# empty log files folder:
rm -r /path/to/project
rm -r /path/to/project
