# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# MEM=55000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R \

library(data.table)
library(R.utils)
library(tidyverse)

rm(list=ls())

# only multithreaded when writing
setDTthreads(4)

path<-"/path/to/project"
dataset<-"macromap"

# provide trait (in this case combination of stimuly + time) as input
args = commandArgs(trailingOnly=TRUE)
conditions<-args[1]

# to test
# conditions<-"CIL_24"

print(conditions)

columpsToKeep<-c("variant","r2","pvalue","molecular_trait_object_id","molecular_trait_id","maf",
                 "gene_id","median_tpm","beta","se","an","ac","chromosome","position","ref","alt","type","rsid")

########################
## NOMINAL DATA:

# get exact file - includes PC number
file<-list.files(paste0(path,dataset,"/raw_data/nominal/",sep=""))
file<-file[grep(conditions,file)]
file<-paste0(path,dataset,"/raw_data/nominal/",file,sep="")

# load file
df1<-fread(file,header=T,sep="\t")


# colnames(df1)<-c("molecular_trait_object_id","Chromosome_phe","TSS","TSS_end","Strand","Total_no_variants_cis","Distance_tss_variant","Variant_in_cis","chromosome","position","position_end","pvalue","beta","Best_hit","se","gene_snp")
# rename some columns
setnames(df1, old = c('Nominal_pvalue','Beta_regression','std.err','Chromosome_var','Pos_variant','Phenotype_ID'), new = c('pvalue','beta','se','chromosome','position','molecular_trait_id'))


# load additional information
snp_info<-fread(paste0(path,dataset,"/raw_data/nominal/snp.info.txt.gz"),head=T)
snp_info$Variant_in_cis<-paste(snp_info$ID,snp_info$CHROM,snp_info$POS,sep="_")
snp_info$variant<-paste(snp_info$CHROM,snp_info$POS,snp_info$REF,snp_info$ALT,sep="_")
snp_info$rsid<-snp_info$variant

snp_info$Alt_freq<-NA
snp_info$Alt_freq[which(snp_info$ALT==snp_info$MINOR)]<-snp_info$MAF[which(snp_info$ALT==snp_info$MINOR)]
snp_info$Alt_freq[which(snp_info$REF==snp_info$MINOR)]<-1-snp_info$MAF[which(snp_info$REF==snp_info$MINOR)]

# four SNP IDs (2 variants) where alleles do not match -but ref = major (minor similar maf)
snp_info$Alt_freq[which(snp_info$ID %in% c("chr6_32663512","chr17_45646564") & is.na(snp_info$Alt_freq))]<-1-snp_info$MAF[which(snp_info$ID %in% c("chr6_32663512","chr17_45646564") & is.na(snp_info$Alt_freq))]

sample<-fread(paste0(path,dataset,"/raw_data/sample_size_macromap.csv"))
sample<-sample[which(sample$condition==conditions),]

snp_info$an<-2*(sample$sample_size)
snp_info$ac<-round((snp_info$an)*(snp_info$Alt_freq))

# rename some columns
setnames(snp_info, old = c('REF','ALT','MAF'), new = c('ref','alt','maf'))

df1<-merge(df1,snp_info[,c("Variant_in_cis","variant","ref","alt","rsid","maf","ac","an")],by="Variant_in_cis")

df1$molecular_trait_object_id<-df1$molecular_trait_id
df1$gene_id<-df1$molecular_trait_id
df1$r2<-NA
df1$median_tpm<-NA
df1$type<-NA

df1<-df1[,..columpsToKeep]
# df1<-df1[!duplicated(df1),]

df1$chromosome<-gsub("chr","",df1$chromosome)
df1$chromosome<-as.numeric(gsub("X","23",df1$chromosome))

df1<-df1[order(chromosome,position,decreasing=F),]

# make new directory if needed:
dir1<-paste(path,dataset,"/nominal/",conditions,"/",sep="")
if (!dir.exists(dir1)) {
  system(paste("mkdir ",dir1,sep=""))
}

# save nominal file, and index:
file_out<-paste(path,dataset,"/nominal/",conditions,"/",conditions,".tsv",sep="")
fwrite(df1,file=file_out,col.names = T,row.names = FALSE,sep="\t",quote = FALSE,na="NA")

# remove existing files
system(paste("rm ",file_out,".gz",sep=""))
system(paste("rm ",file_out,".gz.tbi",sep=""))

system(paste("bgzip ",file_out,sep=""))
system(paste("tabix -f -s13 -b14 -e14 -S1 ",file_out,".gz",sep=""))

# rm(df1)


####################
# create sample size file

# make new directory if needed:
dir1<-paste(path,dataset,"/sample_size_per_tissue/",sep="")
if (!dir.exists(dir1)) {
  system(paste("mkdir ",dir1,sep=""))
}

colnames(sample)<-c("dataset","n")

file.out<-paste(path,dataset,"/sample_size_per_tissue/",sample$dataset,sep="")
fwrite(sample,file.out,col.names=F,row.names=F,quote=F,sep="\t")

########################
# PERMUTED DATA:

file2<-paste0(path,dataset,"/raw_data/permuted/",conditions,".permuted.txt.gz",sep="")
df2<-fread(file2,header=TRUE,sep="\t")

# add new variant name - note some empty rows in df2 (no SNP link to gene) thus merged df2 smaller than one saved by Nikos
df2<-merge(df2,snp_info[,c("ID","variant")],by.x="Best_variant_in_cis",by.y="ID")
rm(snp_info)

# already in the format required by the script
df_final<-data.frame(Phenotype_ID=df2$Phenotype_ID,
                     Chromosome_phe=df2$Chromosome_phe,
                     TSS=df2$TSS,
                     TSS_end=df2$TSS_end,
                     Strand=df2$Strand,
                     Total_no_variants_cis=df2$Total_no_variants_cis,
                     Distance_tss_variant=df2$Distance_tss_variant,
                     Best_variant_in_cis=df2$variant,
                     Chromosome_var=df2$Chromosome_var,
                     Pos_variant=df2$Pos_variant,
                     Pos_variant_end=df2$Pos_variant_end,
                     DF=NA,
                     Dummy=df2$Dummy,
                     Beta_dist_1=df2$Beta_dist_1,
                     Beta_dist_2_number_of_ind_tests=df2$Beta_dist_2_number_of_ind_tests,
                     Nominal_pvalue=df2$Nominal_pvalue,
                     Beta_regression=df2$Beta_regression,
                     Empirical_pvalue=df2$Empirical_pvalue,
                     Corrected_pvalue=df2$Corrected_pvalue,
                     std.err=NA)
rm(df2)

df_final$Chromosome_phe<-gsub("chr","",df_final$Chromosome_phe)
df_final$Chromosome_var<-gsub("chr","",df_final$Chromosome_var)

# make new directory if needed:
dir2<-paste(path,dataset,"/permuted/",sep="")
if (!dir.exists(dir2)) {
  system(paste("mkdir ",dir2,sep=""))
}

# save file, and index:
file_out_2<-paste(path,dataset,"/permuted/",conditions,".permuted.txt",sep="")

fwrite(df_final,file_out_2,col.names = TRUE,row.names = FALSE,sep="\t",quote = FALSE,na="NA")

system(paste("rm ",file_out_2,".gz",sep=""))
system(paste("bgzip ",file_out_2))

q("no")


