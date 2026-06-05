# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# MEM=55000
# bsub -Is -m "modern_hardware" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -n 4 R \

library(data.table)
library(R.utils)
library(dplyr)

# only multithreaded when writing
setDTthreads(4)

dataset<-"hu_2021"

path<-paste0("/path/to/project",dataset,"/")

# provide trait (in this case combination of database + trait) as input
args = commandArgs(trailingOnly=TRUE)
# study<-args[1]
# print(study)
# 
# dataset<-args[2]
# print(dataset)

conditions<-"intestinal_mucosa_colon_ileum"

columpsToKeep<-c("variant","r2","pvalue","molecular_trait_object_id","molecular_trait_id","maf",
                 "gene_id","median_tpm","beta","se","an","ac","chromosome","position","ref","alt","type","rsid")

columpsToKeep_permuted<-c("Phenotype_ID","Chromosome_phe","TSS","TSS_end","Strand","Total_no_variants_cis",
                          "Distance_tss_variant","Best_variant_in_cis","Chromosome_var","Pos_variant","Pos_variant_end","DF",
                          "Dummy","Beta_dist_1","Beta_dist_2_number_of_ind_tests","Nominal_pvalue","Beta_regression",
                          "Empirical_pvalue","Corrected_pvalue","std.err")

## map b37 - b38 file:

df2<-fread(paste0(path,"other/snp_info_b37_b38_tmp.txt.gz"),
           header=TRUE,sep="\t")

# exclude the variants with larger discrepancies with maf in nfe, and those not in gnomad, where the check could not have been carried out
table(df2$col)
# 0        1 
# 4179043    3586

df2<-df2[which(df2$col==0),]
dim(df2)
# [1] 4179043       18

df2$position_b38<-gsub("[0-9]{1,2}:","",df2$variant_b38)
df2$position_b38<-gsub("_.*","",df2$position_b38)

colnames(df2)[2]<-"variant_b37"
df2<-df2[,c("variant_b37","variant_b38","chr","position_b38","Ref","Alt","Alt_freq")]
dim(df2)
df2<-df2[!duplicated(df2),]
dim(df2)


###############################
# original dataset in b37

d<-fread(paste(path,"raw_data/Merged.FDR.txt.gz",sep=""))

d$variant_b37<-paste("chr",d$Chr,"_",d$Pos,"_",d$Allele0,"_",d$Allele1,sep="")

dim(d)
# [1] 44409897       17
dim(d[which(d$rsID %in% df2$variant_b37),])
# [1] 42609976       17

# combine with original dataset:
dim(d)
# [1] 44254957       24

d<-merge(d,df2,by="variant_b37")
dim(d)
# [1] 44254957       24
rm(df2)

###############################################################
# re estimate betas based on alleles in b37 and b38
d$beta_ed<-NA

# no change - beta
d$beta_ed[which(d$Allele1==d$Alt & d$Allele0==d$Ref)]<-d$Beta[which(d$Allele1==d$Alt & d$Allele0==d$Ref)]
dim(d[is.na(d$beta_ed),])
# [1] 11050172    14

# flip - beta*-1
d$beta_ed[which(d$Allele1==d$Ref & d$Allele0==d$Alt)]<- (d$Beta[which(d$Allele1==d$Ref & d$Allele0==d$Alt)])*(-1)
dim(d[is.na(d$beta_ed),])
#[1] 1014   14

# VARIANTS WHERE Ref AND Alt HAVE BEEN SWAPPED -  same beta:

d$beta_ed[which( (d$Allele0=="A") & (d$Ref=="T") & (d$Allele1=="C") & (d$Alt=="G"))]<-d$Beta[which( (d$Allele0=="A") & (d$Ref=="T") & (d$Allele1=="C") & (d$Alt=="G"))]
d$beta_ed[which( (d$Allele0=="A") & (d$Ref=="T") & (d$Allele1=="G") & (d$Alt=="C"))]<-d$Beta[which( (d$Allele0=="A") & (d$Ref=="T") & (d$Allele1=="G") & (d$Alt=="C"))]
d$beta_ed[which( (d$Allele0=="C") & (d$Ref=="G") & (d$Allele1=="A") & (d$Alt=="T"))]<-d$Beta[which( (d$Allele0=="C") & (d$Ref=="G") & (d$Allele1=="A") & (d$Alt=="T"))]
d$beta_ed[which( (d$Allele0=="C") & (d$Ref=="G") & (d$Allele1=="T") & (d$Alt=="A"))]<-d$Beta[which( (d$Allele0=="C") & (d$Ref=="G") & (d$Allele1=="T") & (d$Alt=="A"))]
d$beta_ed[which( (d$Allele0=="G") & (d$Ref=="C") & (d$Allele1=="T") & (d$Alt=="A"))]<-d$Beta[which( (d$Allele0=="G") & (d$Ref=="C") & (d$Allele1=="T") & (d$Alt=="A"))]
d$beta_ed[which( (d$Allele0=="G") & (d$Ref=="C") & (d$Allele1=="A") & (d$Alt=="T"))]<-d$Beta[which( (d$Allele0=="G") & (d$Ref=="C") & (d$Allele1=="A") & (d$Alt=="T"))]
d$beta_ed[which( (d$Allele0=="T") & (d$Ref=="A") & (d$Allele1=="C") & (d$Alt=="G"))]<-d$Beta[which( (d$Allele0=="T") & (d$Ref=="A") & (d$Allele1=="C") & (d$Alt=="G"))]
d$beta_ed[which( (d$Allele0=="T") & (d$Ref=="A") & (d$Allele1=="G") & (d$Alt=="C"))]<-d$Beta[which( (d$Allele0=="T") & (d$Ref=="A") & (d$Allele1=="G") & (d$Alt=="C"))]

dim(d[which(is.na(d$beta_ed)),])
# [1] 10 25


# VARIANTS WHERE Ref AND Alt HAVE BEEN SWAPPED AND FLIPPED - beta*-1:

d$beta_ed[which( (d$Allele0=="A") & (d$Ref=="G") & (d$Allele1=="C") & (d$Alt=="T"))]<-d$Beta[which( (d$Allele0=="A") & (d$Ref=="G") & (d$Allele1=="C") & (d$Alt=="T"))]*(-1)
d$beta_ed[which( (d$Allele0=="A") & (d$Ref=="C") & (d$Allele1=="G") & (d$Alt=="T"))]<-d$Beta[which( (d$Allele0=="A") & (d$Ref=="C") & (d$Allele1=="G") & (d$Alt=="T"))]*(-1)
d$beta_ed[which( (d$Allele0=="C") & (d$Ref=="T") & (d$Allele1=="A") & (d$Alt=="G"))]<-d$Beta[which( (d$Allele0=="C") & (d$Ref=="T") & (d$Allele1=="A") & (d$Alt=="G"))]*(-1)
d$beta_ed[which( (d$Allele0=="C") & (d$Ref=="A") & (d$Allele1=="T") & (d$Alt=="G"))]<-d$Beta[which( (d$Allele0=="C") & (d$Ref=="A") & (d$Allele1=="T") & (d$Alt=="G"))]*(-1)
d$beta_ed[which( (d$Allele0=="G") & (d$Ref=="A") & (d$Allele1=="T") & (d$Alt=="C"))]<-d$Beta[which( (d$Allele0=="G") & (d$Ref=="A") & (d$Allele1=="T") & (d$Alt=="C"))]*(-1)
d$beta_ed[which( (d$Allele0=="G") & (d$Ref=="T") & (d$Allele1=="A") & (d$Alt=="C"))]<-d$Beta[which( (d$Allele0=="G") & (d$Ref=="T") & (d$Allele1=="A") & (d$Alt=="C"))]*(-1)
d$beta_ed[which( (d$Allele0=="T") & (d$Ref=="G") & (d$Allele1=="C") & (d$Alt=="A"))]<-d$Beta[which( (d$Allele0=="T") & (d$Ref=="G") & (d$Allele1=="C") & (d$Alt=="A"))]*(-1)
d$beta_ed[which( (d$Allele0=="T") & (d$Ref=="C") & (d$Allele1=="G") & (d$Alt=="A"))]<-d$Beta[which( (d$Allele0=="T") & (d$Ref=="C") & (d$Allele1=="G") & (d$Alt=="A"))]*(-1)

dim(d[which(is.na(d$beta_ed)),])
# [1]  0 25


# re-estimate ac and an - used later in the functions to re-estimate maf - keep ac so that it refers to alt allele
# 299 samples from 170 patients used in the study:

d$an<-2*(299-(d$MissingSample))
d$ac<-round((d$an)*(d$Alt_freq))

d$r2<-NA
d$pvalue<-d$p_lrt
d$molecular_trait_object_id<-d$ExpressionGene
d$molecular_trait_id<-d$ExpressionGene
d$maf<-pmin(d$Alt_freq,1-d$Alt_freq)
d$gene_id<-d$ExpressionGene
d$median_tpm<-NA
d$beta<-d$beta_ed
d$se<-d$SE
d$chromosome<-d$Chr
d$position<-d$position_b38
d$ref<-d$Ref
d$alt<-d$Alt

d$type<-"indel"
d$type[which(nchar(d$ref)==1 & nchar(d$alt)==1)]<-"SNP"

d$variant_b38<-paste("chr",d$variant_b38,sep="")
d$variant_b38<-gsub(":","_",d$variant_b38)
d$rsid<-d$variant_b38
d$variant<-d$variant_b38


########################
## NOMINAL DATA:

d$position<-as.numeric(d$position)
d$chromosome<-as.numeric(d$chromosome)

d<-d[order(d$chromosome,d$position,decreasing=F),]

d1<-d[,..columpsToKeep]

# make new directory if needed:
dir1<-paste(path,"nominal/",conditions,"/",sep="")
if (!dir.exists(dir1)) {
  system(paste("mkdir ",dir1,sep=""))
}

# save nominal file, and index:
d1$position<-format(d1$position,scientific = FALSE)

file_out<-paste(path,"nominal/",conditions,"/",conditions,".tsv",sep="")
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

# add TSS for each gene - same annotation as in eQTL catalog
gtf<-rtracklayer::import("/path/to/project")
gene<-as.data.frame(gtf)
rm(gtf)

gene$gene_id.2<-gsub("\\.[0-9]{1,2}$","",gene$gene_id)
gene<-gene[which( (gene$gene_id.2 %in% d$molecular_trait_id) & (gene$type=="gene")),c("gene_id.2","start","end","strand")]

gene$TSS<-NA
gene$TSS[which(gene$strand=="+")]<-gene$start[which(gene$strand=="+")]
gene$TSS[which(gene$strand=="-")]<-gene$end[which(gene$strand=="-")]

d<-merge(d,gene[,c("gene_id.2","TSS","strand")],by.x="molecular_trait_id",by.y="gene_id.2",all.x=T)
d$Distance_tss_variant<-as.numeric(d$position)-d$TSS

# get N gene IDs not in the version of gencode used to define TSS - deprecated already in gencode v35
length(gene_ids[which(!gene_ids %in% gene$gene_id.2)])
# [1] 1324
dim(d[which(is.na(d$strand)),])
# [1] 1873434      45
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
df<-d[order(d$FDR,decreasing=F),]
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
                     Nominal_pvalue=df$p_wald,# on which FDR is carried out
                     Beta_regression=df$beta,
                     Empirical_pvalue=NA,
                     Corrected_pvalue=df$FDR,
                     std.err=NA)


df_final<-df_final[which(!is.na(df_final$Corrected_pvalue) & df_final$Corrected_pvalue!=""),]

df_final$Best_variant_in_cis[which(df_final$Best_variant_in_cis==".")] <- as.character(df_final$Dummy[which(df_final$Best_variant_in_cis==".")])
df_final<-df_final[order(df_final$Chromosome_phe,df_final$Pos_variant,decreasing =F),]

print(paste("N FINAL variants in nominal available in permuted: ",nrow(df_final),sep=""))
#[1] "N FINAL variants in nominal available in permuted: 27590"

# make new directory if needed:
dir2<-paste(path,"permuted/",sep="")
if (!dir.exists(dir2)) {
  system(paste("mkdir ",dir2,sep=""))
}

# save file, and index:
file_out_2<-paste(path,"permuted/",conditions,".permuted.txt",sep="")

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

file_out<-paste(path,"variant_info/",conditions,".variant.info.tsv",sep="")

df$position<-format(df$position,scientific = FALSE)
fwrite(df,file=file_out,col.names = F,row.names = FALSE,sep="\t",quote = FALSE)

system(paste("rm ",file_out,".gz",sep=""))
system(paste("bgzip ",file_out,sep=""))
system(paste("tabix -s1 -b2 -e2 ",file_out,".gz",sep=""))

####################
# create sample size file

d<-as.data.frame(matrix(ncol=2))
colnames(d)<-c("dataset","n")

d$dataset<-conditions
d$n<-299

file.out<-paste(path,"sample_size_per_tissue/",d$dataset,sep="")
fwrite(d,file.out,col.names=F,row.names=F,quote=F,sep="\t")

####################

q("no")


