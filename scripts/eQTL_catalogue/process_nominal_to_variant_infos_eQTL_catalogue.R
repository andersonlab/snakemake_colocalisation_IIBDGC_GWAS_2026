# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# MEM=30000
# bsub -Is -m "modern_hardware" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -n 4 R \

# adapted from original file by Nikos in /path/to/project

library(data.table)
library(R.utils)
library(dplyr)

# only multithreaded when writting
setDTthreads(4)

path<-"/path/to/project"

args = commandArgs(trailingOnly=TRUE)
study<-args[1]
dataset<-args[2]

# study<-"QTS000015"
# dataset<-"QTD000116"

columpsToKeep<-c("chromosome","position","rsid","ref","alt","type","ac","an")

file<-paste(path,"nominal/",study,"_",dataset,"/",study,"_",dataset,".tsv.gz",sep="")

df<-fread(file,select=columpsToKeep,header=TRUE)
df$rsid[which(is.na(df$rsid) | (df$rsid==""))]<-paste(df$chromosome,df$position,df$ref,df$alt,sep="_")[which(is.na(df$rsid) | (df$rsid==""))]
df<-df %>% group_by(rsid) %>% distinct()

file_out<-paste(path,"variant_info/",study,"_",dataset,".variant.info.tsv",sep="")

fwrite(df,file=file_out,col.names = F,row.names = FALSE,sep="\t",quote = FALSE)

system(paste("rm ",file_out,".gz",sep=""))
system(paste("bgzip ",file_out,sep=""))
system(paste("tabix -s1 -b2 -e2 ",file_out,".gz",sep=""))

q("no")
