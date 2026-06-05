#!/bin/bash
# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
#BSUB -o sm_logs/snakemake_master-%J-output.log
#BSUB -e sm_logs/snakemake_master-%J-error.log 
#BSUB -q oversubscribed
#BSUB -G your_hpc_group
#BSUB -n 1
#BSUB -M 10000
#BSUB -a "memlimit=True"
#BSUB -R "select[mem>10000] rusage[mem=10000] span[hosts=1]"
#BSUB -J 1

# ${study} defined as an env variable at submission

# Define some params as enviromental variables, done before script execution via:

# export eqtl_study=eQTL_catalogue # TO BE EDITED
# export gwas_study=iibdgc # TO BE EDITED

script_dir=/path/to/project

# worfklow_prefix="coloc_"
group="team152"

workdir=/path/to/project

# Load snakemake and singulatiry
# singularity exec iibdgc_postprocess_10_singularity.sif
# singularity exec iibdgc_postprocess_10_singularity.sif
which singularity

echo ${workdir}results/${eqtl_study}_${gwas_study}/

# Copy config to results
cp ${script_dir}run_coloc_eQTL_per_eqtl_study.smk ${workdir}results/${eqtl_study}_${gwas_study}/

# Copy snakemake cluster config files to the working directory
cp ${script_dir}coloc_config_eQTL_per_eqtl_study.yaml ${workdir}results/${eqtl_study}_${gwas_study}/
cp ${script_dir}cluster_config.yaml ${workdir}results/${eqtl_study}_${gwas_study}/

# Copy main code used in the analyses:
cp ${script_dir}GWAS_run_coloc_eQTL_LF.R ${workdir}results/${eqtl_study}_${gwas_study}/
cp ${script_dir}functions_me_eQTL_bh.R ${workdir}results/${eqtl_study}_${gwas_study}/

# Make new directory for the result files:
mkdir -p ${script_dir}run_coloc_eQTL_per_eqtl_study.smk ${workdir}results/${eqtl_study}_${gwas_study}/coloc_results/

# workdir not required if work launched from right root directory

# Run snakemake
snakemake -j 20000 \
    --latency-wait 90 \
    --use-envmodules \
    --keep-incomplete \
    --keep-going \
    --default-resources threads=1 mem_mb=2000 \
    --directory ${workdir} \
    --cluster-config ${workdir}results/${eqtl_study}_${gwas_study}/coloc_config_eQTL_per_eqtl_study.yaml \
    --cluster-config ${workdir}results/${eqtl_study}_${gwas_study}/cluster_config.yaml \
    --use-singularity \
    --singularity-args "-B /lustre,/software" \
    --keep-going \
    --restart-times 0 \
    --snakefile ${workdir}results/${eqtl_study}_${gwas_study}/run_coloc_eQTL_per_eqtl_study.smk

# bsub -M 2000 -a "memlimit=True" -R "select[mem>2000] rusage[mem=2000] span[hosts=1]" -o sm_logs/snakemake_master-%J-output.log -e sm_logs/snakemake_master-%J-error.log -q oversubscribed -J "snakemake_master_COLOC" < submit_snakemake_BH.sh 
