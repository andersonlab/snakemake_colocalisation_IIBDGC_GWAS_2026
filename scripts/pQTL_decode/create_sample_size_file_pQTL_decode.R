# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#


# # singularity exec iibdgc_postprocess_10_singularity.sif

# MEM=2000
# bsub -Is -M"$MEM" -R"select[model==Intel_Platinum && mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q yesterday R

library(data.table)
library(R.utils)
library(dplyr)
library(rtracklayer)
library(qvalue)

rm(list=ls())

path<-"/path/to/project"


# provide trait (in this case decode or ukb pqtl data)
args = commandArgs(trailingOnly=TRUE)

dataset<-args[1]
print(dataset)

condition<-args[2]
print(condition)

# dataset<-"pQTL_decode"
# condition<-"final_olink_ukb_bi"
# condition<-"final_somascan_smp"

####################
# create sample size file

d3<-as.data.frame(matrix(ncol=2))
colnames(d3)<-c("dataset","n")

d3$dataset<-condition

if (condition=="final_olink_ukb_bi") {
    d3$n<-46218
} else if (condition=="final_somascan_smp") {
    d3$n<-36136
}

# make new directory if needed:
dir4<-paste(path,dataset,"/sample_size_per_condition/",sep="")
if (!dir.exists(dir4)) {
  system(paste("mkdir -p ",dir4,sep=""))
}

file.out<-paste(path,dataset,"/sample_size_per_condition/",d3$dataset,sep="")
fwrite(d3,file.out,col.names=F,row.names=F,quote=F,sep="\t")

####################

q("no")


