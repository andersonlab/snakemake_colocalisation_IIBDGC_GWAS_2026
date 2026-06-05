# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#

# singularity exec iibdgc_postprocess_10_singularity.sif

# # path="/path/to/project"
path="/path/to/project"

# /path/to/project
# /path/to/project

mkdir -p ${path}eQTL_catalogue/nominal/
mkdir -p ${path}eQTL_catalogue/permuted/
mkdir -p ${path}eQTL_catalogue/sample_size_per_condition/
mkdir -p ${path}eQTL_catalogue/log/
mkdir -p ${path}eQTL_catalogue/variant_info/


#######################################################################################################
# 0.- list complete list of studies in eQTL catalog:


# MEM=1000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -q normal R

library(data.table)

path="/path/to/project"

# metadata version 6, 7, and post lustre recovery:
# d<-fread("/path/to/project")
# d<-fread("~/git/snakemake_colocalisation/scripts/eQTL_catalogue/dataset_metadata_eqtl_v7.tsv")

d<-fread("/path/to/project")

names(table(d$study_label))
#  [1] "Alasoo_2018"                "Aygun_2021"                
#  [3] "BLUEPRINT"                  "Bossini-Castillo_2019"     
#  [5] "Braineac2"                  "BrainSeq"                  
#  [7] "CAP"                        "CEDAR"                     
#  [9] "CommonMind"                 "Cytoimmgen"                
# [11] "Fairfax_2012"               "Fairfax_2014"              
# [13] "FUSION"                     "GENCORD"                   
# [15] "GEUVADIS"                   "Gilchrist_2021"            
# [17] "GTEx"                       "HipSci"                    
# [19] "iPSCORE"                    "Jerber_2021"               
# [21] "Kasela_2017"                "Kim-Hellmuth_2017"         
# [23] "Lepik_2017"                 "Naranbhai_2015"            
# [25] "Nathan_2022"                "Nedelec_2016"              
# [27] "OneK1K"                     "OTAR2057_IBDverse"         
# [29] "OTAR2065_neuroinflammation" "OTAR3087_MAGE"             
# [31] "Peng_2018"                  "Perez_2022"                
# [33] "PhLiPS"                     "PISA"                      
# [35] "Quach_2016"                 "Randolph_2021"             
# [37] "ROSMAP"                     "Schmiedel_2018"            
# [39] "Schwartzentruber_2018"      "Steinberg_2020"            
# [41] "Sun_2018"                   "TwinsUK"                   
# [43] "van_de_Bunt_2015"           "Walker_2019"               
# [45] "Young_2019"


names(table(d$tissue_label))
#   [1] "adipose"                                                                   
#   [2] "adipose (visceral)"                                                        
#   [3] "adrenal gland"                                                             
#   [4] "anorectum goblet cell"                                                     
#   [5] "artery (aorta)"                                                            
#   [6] "artery (coronary)"                                                         
#   [7] "artery (tibial)"                                                           
#   [8] "astrocyte"                                                                 
#   [9] "B cell"                                                                    
#  [10] "BEST4+ colonocyte"                                                         
#  [11] "BEST4+ enterocyte"                                                         
#  [12] "blood"                                                                     
#  [13] "brain (amygdala)"                                                          
#  [14] "brain (anterior cingulate cortex)"                                         
#  [15] "brain (caudate)"                                                           
#  [16] "brain (cerebellum)"                                                        
#  [17] "brain (cortex)"                                                            
#  [18] "brain (DLPFC)"                                                             
#  [19] "brain (hippocampus)"                                                       
#  [20] "brain (hypothalamus)"                                                      
#  [21] "brain (nucleus accumbens)"                                                 
#  [22] "brain (putamen)"                                                           
#  [23] "brain (spinal cord)"                                                       
#  [24] "brain (substantia nigra)"                                                  
#  [25] "breast"                                                                    
#  [26] "cartilage"                                                                 
#  [27] "CD16+ monocyte"                                                            
#  [28] "CD4-positive, alpha-beta memory T cell"                                    
#  [29] "CD4-positive, alpha-beta memory T cell, CD45RO-positive"                   
#  [30] "CD4-positive, alpha-beta T cell"                                           
#  [31] "CD4-positive, CD25-positive, alpha-beta regulatory T cell"                 
#  [32] "CD4+ CTL cell"                                                             
#  [33] "CD4+ memory T cell"                                                        
#  [34] "CD4+ T cell"                                                               
#  [35] "CD4+ TCM cell"                                                             
#  [36] "CD4+ TEM cell"                                                             
#  [37] "CD56+ NK cell"                                                             
#  [38] "CD8-positive, alpha-beta memory T cell, CD45RO-positive"                   
#  [39] "CD8+ T cell"                                                               
#  [40] "CD8+ TCM cell"                                                             
#  [41] "CD8+ TEM cell"                                                             
#  [42] "central memory CD4-positive, alpha-beta T cell"                            
#  [43] "classical monocyte"                                                        
#  [44] "colonocyte"                                                                
#  [45] "conventional dendritic cell"                                               
#  [46] "dendritic cell"                                                            
#  [47] "dnT cell"                                                                  
#  [48] "dopaminergic neuron"                                                       
#  [49] "early colonocyte"                                                          
#  [50] "effector memory CD8-positive, alpha-beta T cell"                           
#  [51] "effector memory CD8-positive, alpha-beta T cell, terminally differentiated"
#  [52] "enterochromaffin-like cell"                                                
#  [53] "enterocyte of epithelium proper of ileum"                                  
#  [54] "ependymal cell"                                                            
#  [55] "esophagus (gej)"                                                           
#  [56] "esophagus (mucosa)"                                                        
#  [57] "esophagus (muscularis)"                                                    
#  [58] "fibroblast"                                                                
#  [59] "floor plate progenitor"                                                    
#  [60] "gamma-delta T cell"                                                        
#  [61] "gdT cell"                                                                  
#  [62] "granzyme K-associated CD8 T cell"                                          
#  [63] "heart (atrial appendage)"                                                  
#  [64] "heart (left ventricle)"                                                    
#  [65] "hematopoietic precursor cell"                                              
#  [66] "hepatocyte"                                                                
#  [67] "ileal goblet cell"                                                         
#  [68] "ileum"                                                                     
#  [69] "innate lymphoid cell"                                                      
#  [70] "intermediate monocyte"                                                     
#  [71] "intestinal crypt stem cell of colon"                                       
#  [72] "intestinal crypt stem cell of small intestine"                             
#  [73] "intestinal enteroendocrine cell"                                           
#  [74] "intestinal tuft cell"                                                      
#  [75] "iPSC"                                                                      
#  [76] "kidney (cortex)"                                                           
#  [77] "LCL"                                                                       
#  [78] "liver"                                                                     
#  [79] "lung"                                                                      
#  [80] "macrophage"                                                                
#  [81] "MAIT cell"                                                                 
#  [82] "mast cell"                                                                 
#  [83] "memory B cell"                                                             
#  [84] "microglia"                                                                 
#  [85] "minor salivary gland"                                                      
#  [86] "monocyte"                                                                  
#  [87] "muscle"                                                                    
#  [88] "myofibroblast cell"                                                        
#  [89] "naive B cell"                                                              
#  [90] "naive thymus-derived CD8-positive, alpha-beta T cell"                      
#  [91] "natural killer cell"                                                       
#  [92] "neocortex"                                                                 
#  [93] "neural progenitor"                                                         
#  [94] "neuroblast"                                                                
#  [95] "neuron"                                                                    
#  [96] "neutrophil"                                                                
#  [97] "NK cell"                                                                   
#  [98] "non-classical monocyte"                                                    
#  [99] "ovary"                                                                     
# [100] "pancreas"                                                                  
# [101] "pancreatic islet"                                                          
# [102] "paneth cell of epithelium of small intestine"                              
# [103] "pituitary"                                                                 
# [104] "placenta"                                                                  
# [105] "plasma"                                                                    
# [106] "plasma cell"                                                               
# [107] "plasmablast"                                                               
# [108] "plasmacytoid dendritic cell"                                               
# [109] "platelet"                                                                  
# [110] "prostate"                                                                  
# [111] "rectum"                                                                    
# [112] "sensory neuron"                                                            
# [113] "serotonergic neuron"                                                       
# [114] "sigmoid colon"                                                             
# [115] "skin"                                                                      
# [116] "skin (suprapubic)"                                                         
# [117] "small intestine"                                                           
# [118] "spleen"                                                                    
# [119] "stomach"                                                                   
# [120] "synovium"                                                                  
# [121] "T cell"                                                                    
# [122] "testis"                                                                    
# [123] "Tfh cell"                                                                  
# [124] "Th1 cell"                                                                  
# [125] "Th17 cell"                                                                 
# [126] "Th2 cell"                                                                  
# [127] "thyroid"                                                                   
# [128] "tibial nerve"                                                              
# [129] "tissue-resident macrophage"                                                
# [130] "transverse colon"                                                          
# [131] "Treg memory"                                                               
# [132] "Treg naive"                                                                
# [133] "uterus"                                                                    
# [134] "vagina"



# tissues we want to keep:
tissues<-c("anorectum goblet cell","artery (aorta)","artery (coronary)","artery (tibial)",
"B cell","BEST4+ colonocyte","BEST4+ enterocyte","blood","CD16+ monocyte",
"CD4-positive, alpha-beta memory T cell",
"CD4-positive, alpha-beta memory T cell, CD45RO-positive",
"CD4-positive, alpha-beta T cell",
"CD4-positive, CD25-positive, alpha-beta regulatory T cell",
"CD4+ CTL cell",
"CD4+ memory T cell",
"CD4+ T cell",
"CD4+ TCM cell",
"CD4+ TEM cell",
"CD56+ NK cell",
"CD8-positive, alpha-beta memory T cell, CD45RO-positive",
"CD8+ T cell",
"CD8+ TCM cell",
"CD8+ TEM cell",
"central memory CD4-positive, alpha-beta T cell",
"classical monocyte","colonocyte","dnT cell","early colonocyte",
"effector memory CD8-positive, alpha-beta T cell",
"effector memory CD8-positive, alpha-beta T cell, terminally differentiated",
"enterochromaffin-like cell",
"enterocyte of epithelium proper of ileum",
"ependymal cell",
"esophagus (gej)",
"esophagus (mucosa)",
"esophagus (muscularis)",
"fibroblast",
"gamma-delta T cell",
"gdT cell",
"granzyme K-associated CD8 T cell",
"hematopoietic precursor cell",
"ileal goblet cell",
"ileum","innate lymphoid cell" ,"intermediate monocyte","intestinal crypt stem cell of colon",
"intestinal crypt stem cell of small intestine",
"intestinal enteroendocrine cell",
"intestinal tuft cell",
"iPSC","kidney (cortex)",
"LCL","macrophage","MAIT cell","mast cell","memory B cell",
"monocyte","naive B cell",
"naive thymus-derived CD8-positive, alpha-beta T cell","natural killer cell",
"neutrophil","NK cell","non-classical monocyte","paneth cell of epithelium of small intestine",
"plasma","plasma cell","plasmablast","plasmacytoid dendritic cell","platelet",
"rectum","sigmoid colon","small intestine","stomach","T cell","Tfh cell","Th1 cell",
"Th17 cell","Th2 cell","thyroid","tissue-resident macrophage","transverse colon","Treg memory","Treg naive")

# studies profiling those tissues:
table(d$study_label[which(d$tissue_label %in% tissues)])
    #       Alasoo_2018             BLUEPRINT Bossini-Castillo_2019 
    #                20                    15                     5 
    #               CAP                 CEDAR            Cytoimmgen 
    #                10                     9                    63 
    #      Fairfax_2012          Fairfax_2014               GENCORD 
    #                 1                     4                    15 
    #          GEUVADIS        Gilchrist_2021                  GTEx 
    #                 5                     1                    75 
    #            HipSci               iPSCORE           Jerber_2021 
    #                 5                     5                     5 
    #       Kasela_2017     Kim-Hellmuth_2017            Lepik_2017 
    #                 2                     7                     5 
    #    Naranbhai_2015           Nathan_2022          Nedelec_2016 
    #                 1                    29                    15 
    #            OneK1K     OTAR2057_IBDverse         OTAR3087_MAGE 
    #                23                    70                     5 
    #        Perez_2022                PhLiPS            Quach_2016 
    #                 8                     5                    25 
    #     Randolph_2021        Schmiedel_2018              Sun_2018 
    #                12                    75                     1 
    #           TwinsUK 
    #                10 



# studies profiling those tissues:
names(table(d$study_id[which(d$tissue_label %in% tissues)]))
# 

# expand the analyses to:
# gene counts: Conditional quantile normalisation with cqn using gene length and GC content as covariates followed by inverse normal transformation.
# transcript usage: Transcript usage is calculated by dividing the transcript expression estimates (TPM units) the total expression of all transcripts of the same gene. Transcript usage values (0…1 scale) are further standardised using inverse normal transformation.
# LeafCutter: Normalisation for LeafCutter junction usage values are normalised the same way as txrevise and transcript usage estimates.

table(d$quant_method)
# aptamer       exon         ge leafcutter microarray         tx      txrev 
#       1        109        109        109         18        109        109 



#################################################################
# start with ge and microarray, and expand later to other (which are truncated)
d1<-d[which( (d$tissue_label %in% tissues) & (d$quant_method %in% c("ge","aptamer","microarray"))),]

# exclude study QTS000045-OTAR2057_IBDverse and QTS000044-OTAR3087_MAGE - not downloadable yet:
d1<-d1[which(!d1$study_id %in% c("QTS000045","QTS000044")),]

dim(d1)
# [1] 224  11

table(d1$study_label)
    #       Alasoo_2018             BLUEPRINT Bossini-Castillo_2019 
    #                 4                     3                     1 
    #               CAP                 CEDAR            Cytoimmgen 
    #                 2                     9                    63 
    #      Fairfax_2012          Fairfax_2014               GENCORD 
    #                 1                     4                     3 
    #          GEUVADIS        Gilchrist_2021                  GTEx 
    #                 1                     1                    15 
    #            HipSci               iPSCORE           Jerber_2021 
    #                 1                     1                     5 
    #       Kasela_2017     Kim-Hellmuth_2017            Lepik_2017 
    #                 2                     7                     1 
    #    Naranbhai_2015           Nathan_2022          Nedelec_2016 
    #                 1                    29                     3 
    #            OneK1K            Perez_2022                PhLiPS 
    #                23                     8                     1 
    #        Quach_2016         Randolph_2021        Schmiedel_2018 
    #                 5                    12                    15 
    #          Sun_2018               TwinsUK 
    #                 1                     2 


d1$study_id
d1$dataset_id

write.table(d1$study_id,paste0(path,"eQTL_catalogue/list_study_ids"),col.names=F,row.names=F,quote=F)
write.table(d1$dataset_id,paste0(path,"eQTL_catalogue/list_dataset_ids"),col.names=F,row.names=F,quote=F)

q("no")

#######################################################################################################
# 1.- reformat nominal summary stats file, and save in new location preserving the format Nikos uses:
# recreate path hierarchies by Nikos

path="/path/to/project"

studies=($(cat ${path}/eQTL_catalogue/list_study_ids))
dataset=($(cat ${path}/eQTL_catalogue/list_dataset_ids))

# length array
echo ${#dataset[@]}
# 224

for i in {0..223}
do
echo ${studies[${i}]}_${dataset[${i}]} && ls -lh /path/to/project
done

#######################################################

MEM=90000 # most traits could run with 30000, others (ge) up to 90000

for i in {0..223}
do
bsub -J"eQTL" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path}eQTL_catalogue/log/${studies[${i}]}_${dataset[${i}]}_process_stderr \
-o ${path}eQTL_catalogue/log/${studies[${i}]}_${dataset[${i}]}_process_stdout \
"Rscript ~/git/snakemake_colocalisation/scripts/eQTL_catalogue/process_eQTL_catalogue_nominal_and_permuted_summary_stats.R ${studies[${i}]} ${dataset[${i}]} > \
${path}eQTL_catalogue/log/process_eQTL_catalogue_nominal_and_permuted_summary_stats_${studies[${i}]}_${dataset[${i}]}.Rout"
done


for i in {0..223}
do
echo ${studies[${i}]}_${dataset[${i}]} && echo ${i} && tail -50 ${path}eQTL_catalogue/log/${studies[${i}]}_${dataset[${i}]}_process_stdout | grep -E "Successfully|Exited"
done

for i in {0..223}
do
echo ${studies[${i}]}_${dataset[${i}]} && tail -50 ${path}eQTL_catalogue/log/process_eQTL_catalogue_nominal_and_permuted_summary_stats_${studies[${i}]}_${dataset[${i}]}.Rout
done


##############################################
# 2.- Create sample size per trait files:

MEM=100

for i in {0..223}
do
bsub -J"eQTL2" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q "small" \
-e ${path}eQTL_catalogue/log/${studies[${i}]}_${dataset[${i}]}_sample_size_stderr \
-o ${path}eQTL_catalogue/log/${studies[${i}]}_${dataset[${i}]}_sample_size_stdout \
"Rscript ~/git/snakemake_colocalisation/scripts/eQTL_catalogue/create_sample_size_file_eQTL_catalogue.R ${studies[${i}]} ${dataset[${i}]} > \
${path}eQTL_catalogue/log/create_sample_size_file_eQTL_catalogue_${studies[${i}]}_${dataset[${i}]}.Rout"
done

for i in {0..223}
do
echo ${studies[${i}]}_${dataset[${i}]} && tail -50  ${path}eQTL_catalogue/log/${studies[${i}]}_${dataset[${i}]}_sample_size_stdout | grep -E "Successfully|Exited"
done


##############################################
# 3.- generate variant info files - requires step 1 to be completed

MEM=30000 

for i in {0..223}
do
bsub -J"eQTL3" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q "normal" \
-e ${path}eQTL_catalogue/log/${studies[${i}]}_${dataset[${i}]}_variant_info_stderr \
-o ${path}eQTL_catalogue/log/${studies[${i}]}_${dataset[${i}]}_variant_info_stdout \
"Rscript  ~/git/snakemake_colocalisation/scripts/eQTL_catalogue/process_nominal_to_variant_infos_eQTL_catalogue.R ${studies[${i}]} ${dataset[${i}]} > \
${path}eQTL_catalogue/log/process_nominal_to_variant_infos_eQTL_catalogue_${studies[${i}]}_${dataset[${i}]}.Rout"
done


for i in {0..223}
do
echo ${studies[${i}]}_${dataset[${i}]} && tail -50 ${path}eQTL_catalogue/log/${studies[${i}]}_${dataset[${i}]}_variant_info_stdout | grep "Successfully completed."
done


##############################################
# 5.- process the permuted files - already done in step 1


##############################################
# 6.-  final check that all files are in place:


for i in {0..223}
do
echo ${studies[${i}]}_${dataset[${i}]} && \
ls -la ${path}eQTL_catalogue/nominal/${studies[${i}]}_${dataset[${i}]}/${studies[${i}]}_${dataset[${i}]}.tsv.gz && \
ls -la ${path}eQTL_catalogue/nominal/${studies[${i}]}_${dataset[${i}]}/${studies[${i}]}_${dataset[${i}]}.tsv.gz.tbi && \
ls -la ${path}eQTL_catalogue/permuted/${studies[${i}]}_${dataset[${i}]}.permuted.txt.gz && \
ls -la ${path}eQTL_catalogue/sample_size_per_condition/${studies[${i}]}_${dataset[${i}]} && \
ls -la ${path}eQTL_catalogue/variant_info/${studies[${i}]}_${dataset[${i}]}.variant.info.tsv.gz && \
ls -la ${path}eQTL_catalogue/variant_info/${studies[${i}]}_${dataset[${i}]}.variant.info.tsv.gz.tbi
done


##############################################
# 6.-  create a file listing the conditions:

for i in {0..223}
do
echo ${studies[${i}]}_${dataset[${i}]} >>  ${path}eQTL_catalogue/list_conditions.txt
done

##############################################
# 7.-  submit coloc jobs 

export eqtl_study=eQTL_catalogue
export gwas_study=iibdgc

# TO BE EDITED"
# /path/to/project 
# eqtl_study variable NEEDS to be set up

mkdir -p /path/to/project

bsub -M 2000 -a "memlimit=True" -R "select[mem>2000] rusage[mem=2000] span[hosts=1]" \
-o /path/to/project \
-e /path/to/project \
-q oversubscribed -J "snakemake_master_COLOC" < /path/to/project

# Job <270301> is submitted to queue <oversubscribed>.

logs/cluster/run_coloc/output_dir=/path/to/project
# evaluate if all jobs are completed correctly:
tail -50 /path/to/project | grep "Successfully"
head -10  /path/to/project
# job          count
# ---------  -------
# run_all          1
# run_coloc    13440

# number much lower in 2nd attempts - to resumbit larger n chunk with no genes, to get full coloc.DONE 
# Job stats:
# job          count
# ---------  -------
# run_all          1
# run_coloc     1231
# total         1232

# number much lower in 2nd 3rd attempts - to resumbit remove multiallelic sites from GWAS 
# Building DAG of jobs...
# Using shell: /usr/local/bin/bash
# Provided cluster nodes: 20000
# Job stats:
# job          count
# ---------  -------
# run_all          1
# run_coloc      164
# total          165

pheno=(IBD CD UC)

for ph in ${pheno[@]}
do
echo ${ph} && ls -la /path/to/project | wc -l 
done
# IBD
# 4479
# CD
# 4479
# UC
# 4479

# evaluate integrity of files:
pheno=(IBD CD UC)


rm /path/to/project

for ph in ${pheno[@]}
do
echo ${ph} >> /path/to/project && \
for i in {0..223}
do 
for chunck in {1..20}
do
echo ${studies[${i}]}_${dataset[${i}]} >> /path/to/project && \
head -1 /path/to/project | wc -w >> /path/to/project
tail -1 /path/to/project | wc -w >> /path/to/project
done
done
done

less /path/to/project

# combine all data:

export eqtl_study=eQTL_catalogue
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
