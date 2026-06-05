# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
## create file with sample size per tissue, like existing ones for GTEx
# lf9@hpc-login:/path/to/project more Adipose_Subcutaneous 
# Adipose_Subcutaneous	581
# MEM=2000
# bsub -Is -m "modern_hardware" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -n 4 R \


library(data.table)
library(R.utils)

# only multithreaded when writting
setDTthreads(4)

path<-"/path/to/project"

args = commandArgs(trailingOnly=TRUE)
study<-args[1]
dataset<-args[2]


# dataset<-gsub("_ge_.*","",trait)
# dataset<-gsub("_microarray_.*","",dataset)


d<-fread("/path/to/project")
d<-d[which(d$study_id==study & d$dataset_id==dataset),]

df<-as.data.frame(matrix(ncol=2,nrow=1))
colnames(df)<-c("dataset","n")
df$dataset<-paste(study,dataset,sep="_")
df$n<-d$sample_size

file.out<-paste(path,"eQTL_catalogue/sample_size_per_tissue/",df$dataset,sep="")

fwrite(df,file.out,col.names=F,row.names=F,quote=F,sep="\t")

q("no")
