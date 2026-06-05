# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# MEM=30000
# bsub -Is -M"$MEM" -R"select[model==Intel_Platinum && mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group R

# adapted from original file by Nikos in /path/to/project

library(data.table)
library(R.utils)
library(dplyr)

rm(list=ls())

# only multithreaded when writting
setDTthreads(4)

path<-"/path/to/project"
dataset<-"macromap"

args = commandArgs(trailingOnly=TRUE)
conditions<-args[1]

columnsToKeep<-c("chromosome","position","rsid","ref","alt","type","ac","an")

df<-fread(paste0(path,dataset,"/raw_data/nominal/snp.info.txt.gz"),head=T)
setnames(df, old = c('REF','ALT','MAF',"POS"), new = c('ref','alt','maf','position'))

df$chromosome<-gsub("chr","",df$CHROM)
df$rsid<-paste(df$CHROM,df$position,df$ref,df$alt,sep="_")

df$type<-"indel"
df$type[which(nchar(df$ref)==1 & nchar(df$alt)==1)]<-"SNP"

# add sample size:
d<-fread(paste0(path,dataset,"/raw_data/sample_size_macromap.csv"))
colnames(d)<-c("dataset","n")
d<-as.data.frame(d)
d<-d[which(d$dataset==conditions),]

# re-estimate ac and an - used later in the functions to re-estimate maf:
df$an<-d$n*2
df$ac<-round((df$an)*(df$maf))

df<-df[,..columnsToKeep]

df<-df %>% group_by(rsid) %>% distinct()

file_out<-paste(path,dataset,"/variant_info/",conditions,".variant.info.tsv",sep="")

fwrite(df,file=file_out,col.names = F,row.names = FALSE,sep="\t",quote = FALSE)

system(paste("rm ",file_out,".gz",sep=""))
system(paste("bgzip ",file_out,sep=""))
system(paste("tabix -s1 -b2 -e2 ",file_out,".gz",sep=""))

q("no")
