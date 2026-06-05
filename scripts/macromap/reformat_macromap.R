# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
## original data from macromap stored in Gaffney's folder - but HGI gave us granted access: 

# from Allan Daly:
# 
# I sorted the access to the other file
# /path/to/project
# I did this by changing the permissions to the gaffney team folder to allow traversing but still not allowing reading or writing to the folders 
# in the path unless the permissions of the files/folders allows. eg. you still cannot list /path/to/project . 
# I assume that should be fine but it does rely on the people in team gaffney setting the permissions correctly within the team folder rather than 
# just relying on the top folder permissions to block access.

# retrieve data from Zenodo:

# singularity exec iibdgc_postprocess_10_singularity.sif


path="/path/to/project"
dataset="macromap"

mkdir -p ${path}${dataset}/raw_data/
mkdir -p ${path}${dataset}/other/
mkdir -p ${path}${dataset}/log/
mkdir -p ${path}${dataset}/nominal/
mkdir -p ${path}${dataset}/permuted/
mkdir -p ${path}${dataset}/sample_size_per_condition/
mkdir -p ${path}${dataset}/variant_info/

cd ${path}${dataset}/raw_data/

# get data from zenodo:

mkdir permuted/
tar -xvf MacroMap_eQTL_permuted_summary_stats.tar -C permuted/

mkdir nominal/
tar -xvf MacroMap_eQTL_nominal_summary_stats.tar -C nominal/

# sample size file provided by Omar (supplementary table from his paper)


#######################################################################################################
# 1.- reformat nominal summary stats file, and save in new location preserving the format Nikos uses:
# recreate path hierarchies by Nikos


less ${path}${dataset}/raw_data/nominal/README
# # eQTL Summary Statistics
# 
# This repository contains nominal summary statistics of MacroMap eQTL analysis. 
# The summary statistics represent associations between genetic variants (SNPs) and gene expression levels across different conditions. 
# Each file in the repository represents a specific condition at a specific time point. 
# 
# ## Files
# 
# The repository includes the following files:
#   
# CIL_24_1MB_PC40_all.stderr.txt.gz
# CIL_6_1MB_PC40_all.stderr.txt.gz
# Ctrl_24_1MB_PC40_all.stderr.txt.gz
# Ctrl_6_1MB_PC40_all.stderr.txt.gz
# IFNB_24_1MB_PC40_all.stderr.txt.gz
# IFNB_6_1MB_PC35_all.stderr.txt.gz
# IFNG_24_1MB_PC35_all.stderr.txt.gz
# IFNG_6_1MB_PC35_all.stderr.txt.gz
# IL4_24_1MB_PC40_all.stderr.txt.gz
# IL4_6_1MB_PC50_all.stderr.txt.gz
# LIL10_24_1MB_PC35_all.stderr.txt.gz
# LIL10_6_1MB_PC50_all.stderr.txt.gz
# MBP_24_1MB_PC40_all.stderr.txt.gz
# MBP_6_1MB_PC50_all.stderr.txt.gz
# P3C_24_1MB_PC40_all.stderr.txt.gz
# P3C_6_1MB_PC50_all.stderr.txt.gz
# PIC_24_1MB_PC40_all.stderr.txt.gz
# PIC_6_1MB_PC50_all.stderr.txt.gz
# Prec_D0_1MB_PC35_all.stderr.txt.gz
# Prec_D2_1MB_PC35_all.stderr.txt.gz
# R848_24_1MB_PC40_all.stderr.txt.gz
# R848_6_1MB_PC50_all.stderr.txt.gz
# sLPS_24_1MB_PC40_all.stderr.txt.gz
# sLPS_6_1MB_PC50_all.stderr.txt.gz  
# snp.info.txt.gz
# 
# ## SNP Information
# The file snp.info.txt.gz contains additional information about the SNPs analyzed in the eQTL analysis. 
# This file provides details such as major and minor alleles, as well as the reference (REF) and alternate (ALT) alleles for each SNP. 
# The beta values in the summary statistics are calculated based on the REF/ALT alleles.
# 
# Please note that these summary statistics are nominal and not adjusted for multiple testing. 
# For any questions or inquiries, please contact author_contact


#########################################
# after talking with Nikos, headers as in:

less /path/to/project 
# /path/to/project
# /path/to/project
# /path/to/project
# /path/to/project
# /path/to/project
# /path/to/project
# /path/to/project
# /path/to/project
# /path/to/project
# /path/to/project
# /path/to/project
# /path/to/project


studies=(IFNB_6 IFNB_24 P3C_6 P3C_24 CIL_6 CIL_24 IFNG_6 IFNG_24 PIC_6 PIC_24 IL4_6 IL4_24 Prec_D0 Prec_D2 Ctrl_6 Ctrl_24 LIL10_6 LIL10_24 R848_6 R848_24 MBP_6 MBP_24 sLPS_6 sLPS_24)

for s in ${studies[@]}
do 
echo $s && \
ls -la  ${path}${dataset}/raw_data/permuted/${s}.permuted.txt.gz && \
ls -la  ${path}${dataset}/raw_data/nominal/${s}_1MB_*_all.stderr.txt.gz
done


#######################################################################################################
# 1.- reformat nominal summary stats file, and save in new location preserving the format Nikos uses:
# recreate path hierarchies by Nikos

MEM=55000
  
for i in ${studies[@]}
do
bsub -J"eQTL" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path}macromap/log/${i}_process_stderr \
-o ${path}macromap/log/${i}_process_stdout \
"Rscript ~/git/snakemake_colocalisation/scripts/macromap/process_eqtl_macromap_nominal_and_permuted_summary_stats.R ${i} > \
${path}macromap/log/process_eqtl_macromap_nominal_and_permuted_summary_stats_${i}.Rout"
done

# submitted 23 June


for i in ${studies[@]}
do
echo ${i} && tail -50 ${path}macromap/log/${i}_process_stdout | grep -E "Successfully|Exited"
done

##############################################
# 2.- Create sample size per trait files - already done in step 1


##############################################
# 3.- generate variant info files - requires step 1 to be completed

MEM=15000 

for i in ${studies[@]}
do
bsub -J"eQTL" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path}macromap/log/${i}_info_process_stderr \
-o ${path}macromap/log/${i}_info_process_stdout \
"Rscript ~/git/snakemake_colocalisation/scripts/macromap/process_nominal_to_variant_infos_macromap.R ${i} > \
${path}macromap/log/process_nominal_to_variant_infos_macromap_${i}.Rout"
done

# submitted

for i in ${studies[@]}
do
echo ${i} && tail -50 ${path}macromap/log/${i}_info_process_stdout | grep -E "Successfully|Exited"
done

##############################################
# 4.- process the permuted files - already done in step 1


##############################################
# 5.-  final check that all files are in place:


for i in ${studies[@]}
do
echo ${i} && \
ls -la ${path}macromap/nominal/${i}/${i}.tsv.gz && \
ls -la ${path}macromap/nominal/${i}/${i}.tsv.gz.tbi && \
ls -la ${path}macromap/permuted/${i}.permuted.txt.gz && \
ls -la ${path}macromap/sample_size_per_condition/${i} && \
ls -la ${path}macromap/variant_info/${i}.variant.info.tsv.gz && \
ls -la ${path}macromap/variant_info/${i}.variant.info.tsv.gz.tbi
done

# COMPLETED

##############################################
# 6.-  create a file listing the conditions:

studies=(IFNB_6 IFNB_24 P3C_6 P3C_24 CIL_6 CIL_24 IFNG_6 IFNG_24 PIC_6 PIC_24 IL4_6 IL4_24 Prec_D0 Prec_D2 Ctrl_6 Ctrl_24 LIL10_6 LIL10_24 R848_6 R848_24 MBP_6 MBP_24 sLPS_6 sLPS_24)

for i in ${studies[@]}
do 
echo ${i} >>  ${path}macromap/list_conditions.txt
done

##############################################
# 7.-  submit coloc jobs 

# edit /path/to/project to incidate eqtl study

eqtl_study=macromap
gwas_study=iibdgc

mkdir -p /path/to/project

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
# run_coloc     1440

ls -la /path/to/project | wc -l 
# 1440

# evaluate integrity of files:

studies=(IFNB_6 IFNB_24 P3C_6 P3C_24 CIL_6 CIL_24 IFNG_6 IFNG_24 PIC_6 PIC_24 IL4_6 IL4_24 Prec_D0 Prec_D2 Ctrl_6 Ctrl_24 LIL10_6 LIL10_24 R848_6 R848_24 MBP_6 MBP_24 sLPS_6 sLPS_24)
pheno=(IBD CD UC)

for ph in ${pheno[@]}
do
echo ${ph} >> /path/to/project && for s in ${studies[@]}
do 
for chunck in {1..20}
do
echo ${s} >> /path/to/project && \
head -1 /path/to/project | wc -w >> /path/to/project
tail -1 /path/to/project | wc -w >> /path/to/project
done
done
done

less /path/to/project

# combine all data:

export eqtl_study=macromap
export gwas_study=iibdgc

path="/path/to/project"
MEM=200

bsub -J"eQTL" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path}macromap/log/join_coloc_output.stderr \
-o ${path}macromap/log/join_coloc_output.stdout \
"bash /path/to/project"


unset eqtl_study
unset gwas_study

