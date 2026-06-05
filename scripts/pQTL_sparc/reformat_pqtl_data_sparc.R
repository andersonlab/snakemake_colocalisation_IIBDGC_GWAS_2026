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

dataset=pQTL_sparc

mkdir -p ${path}${dataset}/nominal/
mkdir -p ${path}${dataset}/permuted/
mkdir -p ${path}${dataset}/sample_size_per_condition/
mkdir -p ${path}${dataset}/log/
mkdir -p ${path}${dataset}/variant_info/


#######################################################################################################
# 0.- list raw data provided by Kyle:

cd /path/to/project

ls pqtl_results_chr1_* | wc -l
# 2921

ls pqtl_results_chr1_* | sed 's/pqtl_results_chr1_batch[0-9]*_//g' | sed 's/.regenie//g' > \
/path/to/project

# move each protein to one directory:

# exclude encode, remaining encode, Soskic et al, Zhang et al, and immune cell atlas, exclude as well vahedi, not part of baseline
proteins=($(cat /path/to/project 

for pt in ${proteins[@]}
do echo ${pt}
done

echo ${#proteins[@]}
# 2921

cd /path/to/project

for pt in ${proteins[@]}
do 
echo ${pt} && ls -la pqtl_results_chr*_batch*_${pt}.regenie | wc -l
done


for pt in ${proteins[@]}
do 
mkdir ${pt} && \
mv pqtl_results_chr*_batch*_${pt}.regenie ${pt}/
done
done

#######################################################################################################
# 1.- reformat nominal summary stats file, and save in new location preserving the format Nikos uses:
# recreate path hierarchies by Nikos

path="/path/to/project"

MEM=35000 # most traits could run with 30000, others (ge) up to 90000


dataset=pQTL_sparc
conditions=(plasma_ibd_patients)

for condition in ${conditions[@]}
do
bsub -J"eQTL" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q long \
-e ${path}pQTL_decode/log/${dataset}_${condition}_process_stderr \
-o ${path}pQTL_decode/log/${dataset}_${condition}_process_stdout \
"Rscript ~/git/snakemake_colocalisation/scripts/pQTL_sparc/process_pQTL_sparc_nominal_and_permuted_summary_stats.R ${dataset} ${condition} > \
${path}pQTL_sparc/log/process_pQTL_sparc_nominal_and_permuted_summary_stats_${dataset}_${condition}.Rout"
done

# CONTINUE HERE!!!!



for condition in ${conditions[@]}
do
echo ${condition} && echo ${i} && tail -50  ${path}pQTL_decode/log/${dataset}_${condition}_process_stdout | grep -E "Successfully|Exited"
done

for condition in ${conditions[@]}
do
echo ${condition} && tail -50 ${path}${dataset}/log/process_${dataset}_nominal_and_permuted_summary_stats_${condition}.Rout
done


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

bsub -M 2000 -a "memlimit=True" -R "select[mem>2000] rusage[mem=2000] span[hosts=1]" \
-o /path/to/project \
-e /path/to/project \
-q oversubscribed -J "snakemake_master_COLOC" < /path/to/project

# Job <<572940>> is submitted to queue <oversubscribed>.

# evaluate if all jobs are completed correctly:
tail -50 /path/to/project | grep "Successfully"
head -10  /path/to/project
# job          count
# ---------  -------
# run_all          1
# run_coloc       60


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
