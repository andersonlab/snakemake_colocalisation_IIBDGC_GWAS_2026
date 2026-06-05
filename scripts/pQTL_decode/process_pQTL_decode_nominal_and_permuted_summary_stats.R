# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#


# # singularity exec iibdgc_postprocess_10_singularity.sif

# MEM=80000
# bsub -Is -M"$MEM" -R"select[model==Intel_Platinum && mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q yesterday R

library(data.table)
library(R.utils)
library(dplyr)
library(rtracklayer)
library(qvalue)

rm(list=ls())
path_proteomics<-"/path/to/project"
path_gwas<-"/path/to/project"
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

# retrieve cis data for all genes:
d<-fread(paste0(path,dataset,"/tmp/",condition,".txt.gz"))

# exclude any duplicated variant_gene pair:
d<-d[!duplicated(d$variant_gene),]

columpsToKeep<-c("variant","r2","pvalue","molecular_trait_object_id","molecular_trait_id","maf",
                 "gene_id","median_tpm","beta","se","an","ac","chromosome","position","ref","alt","type","rsid")

columpsToKeep_permuted<-c("Phenotype_ID","Chromosome_phe","TSS","TSS_end","Strand","Total_no_variants_cis",
                          "Distance_tss_variant","Best_variant_in_cis","Chromosome_var","Pos_variant","Pos_variant_end","DF",
                          "Dummy","Beta_dist_1","Beta_dist_2_number_of_ind_tests","Nominal_pvalue","Beta_regression",
                          "Empirical_pvalue","Corrected_pvalue","std.err")


###############################


########################
## NOMINAL DATA:

d$position<-as.numeric(d$position)
d$chromosome<-as.numeric(gsub("chr","",d$Chrom))


d<-d[order(d$chromosome,d$position,decreasing=F),]

d1<-d[,..columpsToKeep]

# make new directory if needed:
dir1<-paste(path,dataset,"/nominal/",condition,"/",sep="")
if (!dir.exists(dir1)) {
  system(paste("mkdir -p ",dir1,sep=""))
}

# save nominal file, and index:
d1$position<-format(d1$position,scientific = FALSE)

file_out<-paste(path,dataset,"/nominal/",condition,"/",condition,".tsv",sep="")
fwrite(d1,file=file_out,col.names = T,row.names = FALSE,sep="\t",quote = FALSE,na="NA")

# remove existing files
system(paste("rm ",file_out,".gz",sep=""))
system(paste("rm ",file_out,".gz.tbi",sep=""))

system(paste("bgzip ",file_out,sep=""))
system(paste("tabix -f -s13 -b14 -e14 -S1 ",file_out,".gz",sep=""))

rm(d1)

#############################
# permuted file, for each gene id keep the most significant SNP:

gene_ids<-names(table(d$molecular_trait_id))
length(gene_ids)
# [1] 27590

d$Distance_tss_variant<-as.numeric(d$position)-d$TSS

# # get N gene IDs not in the version of gencode used to define TSS - deprecated already in gencode v35
# length(gene_ids[which(!gene_ids %in% gene$gene_id.2)])
# # [1] 1324
# dim(d[which(is.na(d$strand)),])
# # [1] 1873434      45
# d<-d[which(!is.na(d$strand)),]

# create field d$n_variants - number of variants in cis:

gene_ids<-names(table(d$molecular_trait_id))
length(gene_ids)
#[1] 27590

variants_cis<-as.data.frame(table(d$molecular_trait_id))
colnames(variants_cis)<-c("molecular_trait_id","n_variants")

dim(d)
d<-merge(d,variants_cis,by="molecular_trait_id")
dim(d)

# keep only the most significant variant:
df<-d[order(d$qval,decreasing=F),]
df<-df[!duplicated(df$molecular_trait_id),]


# combine df and permuted data to replicate Nikos file, see how this file is used downstream by coloc pipeline
df_final<-data.frame(Phenotype_ID=df$molecular_trait_id,
                     Chromosome_phe=df$chromosome,
                     TSS=df$TSS,
                     TSS_end=df$TSS,
                     Strand=df$strand,
                     Total_no_variants_cis=df$n_variants,
                     Distance_tss_variant=df$Distance_tss_variant,
                     Best_variant_in_cis=df$variant,
                     Chromosome_var=df$chromosome,
                     Pos_variant=df$position,
                     Pos_variant_end=df$position,
                     DF=NA,
                     Dummy=df$variant,
                     Beta_dist_1=NA,
                     Beta_dist_2_number_of_ind_tests=NA,
                     Nominal_pvalue=df$pvalue,# on which FDR is carried out
                     Beta_regression=df$beta,
                     Empirical_pvalue=NA,
                     Corrected_pvalue=df$qval,
                     std.err=NA)


df_final<-df_final[which(!is.na(df_final$Corrected_pvalue) & df_final$Corrected_pvalue!=""),]

df_final$Best_variant_in_cis[which(df_final$Best_variant_in_cis==".")] <- as.character(df_final$Dummy[which(df_final$Best_variant_in_cis==".")])
df_final<-df_final[order(df_final$Chromosome_phe,df_final$Pos_variant,decreasing =F),]

print(paste("N FINAL variants in nominal available in permuted: ",nrow(df_final),sep=""))
#[1] "N FINAL variants in nominal available in permuted: 27590"

# make new directory if needed:
dir2<-paste(path,dataset,"/permuted/",sep="")
if (!dir.exists(dir2)) {
  system(paste("mkdir ",dir2,sep=""))
}

# save file, and index:
file_out_2<-paste(path,dataset,"/permuted/",condition,".permuted.txt",sep="")

fwrite(df_final,file_out_2,col.names = TRUE,row.names = FALSE,sep="\t",quote = FALSE,na="NA")

system(paste("rm ",file_out_2,".gz",sep=""))
system(paste("bgzip ",file_out_2))
rm(df_final,df)

####################
# create variant info files:

columnsToKeep<-c("chromosome","position","rsid","ref","alt","type","ac","an")

df<-d[,..columnsToKeep]
df<-df[order(df$chromosome,df$position,decreasing=F),]

df<-df %>% group_by(rsid) %>% distinct()

df<-df[!duplicated(df$rsid),]

# make new directory if needed:
dir3<-paste(path,dataset,"/variant_info/",sep="")
if (!dir.exists(dir3)) {
  system(paste("mkdir ",dir3,sep=""))
}

file_out<-paste(path,dataset,"/variant_info/",condition,".variant.info.tsv",sep="")

df$position<-format(df$position,scientific = FALSE)
fwrite(df,file=file_out,col.names = F,row.names = FALSE,sep="\t",quote = FALSE)

system(paste("rm ",file_out,".gz",sep=""))
system(paste("bgzip ",file_out,sep=""))
system(paste("tabix -s1 -b2 -e2 ",file_out,".gz",sep=""))

####################
# create sample size file

d3<-as.data.frame(matrix(ncol=2))
colnames(d3)<-c("dataset","n")

d3$dataset<-condition
d3$n<-max(d$N)

# make new directory if needed:
dir4<-paste(path,dataset,"/sample_size_per_condition/",sep="")
if (!dir.exists(dir4)) {
  system(paste("mkdir -p ",dir4,sep=""))
}

file.out<-paste(path,dataset,"/sample_size_per_condition/",d3$dataset,sep="")
fwrite(d3,file.out,col.names=F,row.names=F,quote=F,sep="\t")

####################

q("no")


