# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# # singularity exec iibdgc_postprocess_10_singularity.sif
# MEM=45000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -n 2 R \

library(data.table)
library(R.utils)

rm(list=ls())

# only multithreaded when writing
setDTthreads(4)

path<-"/path/to/project"

# provide trait (in this case combination of database + trait) as input
args = commandArgs(trailingOnly=TRUE)
study<-args[1]
print(study)

dataset<-args[2]
print(dataset)

# study<-"QTS000036"
# dataset<-"QTD000585"

columpsToKeep<-c("variant","r2","pvalue","molecular_trait_object_id","molecular_trait_id","maf",
                 "gene_id","median_tpm","beta","se","an","ac","chromosome","position","ref","alt","type","rsid")


d<-fread("/path/to/project")
d<-d[which(d$study_id==study & d$dataset_id==dataset),]

conditions<-paste(d$study_id,d$dataset_id,sep="_")
print(conditions)

########################
## NOMINAL DATA:

if (d$quant_method=="ge") {
  file<-paste("/path/to/project",study,"/",dataset,"/",dataset,".all.tsv.gz",sep="")
} else {
  file<-paste("/path/to/project",study,"/",dataset,"/",dataset,".cc.tsv.gz",sep="")
} 

df1<-fread(file,header=TRUE,sep="\t")

# some varinats p-value is already 0 in original file, causing an issue later - modify to min value in dataset other than 0:
# min_pval<-min(min(df1$pvalue[which(df1$pvalue!=0)]),1E-322)
# df1$pvalue[which(df1$pvalue==0)]<-min_pval

# some duplicated variants (diff variants ID) with same rsid (expected, multiallelic ones), and viceversa (?) and variant ids - otherwise the same - rename:
df1$rsid<-df1$variant
print(paste("N tests - with potential duplicates:",nrow(df1)))

df1<-df1[!duplicated(df1),]
print(paste("N tests - without potential duplicates:",nrow(df1)))


########################
# PERMUTED DATA:

file2<-paste("/path/to/project",study,"/",dataset,"/",dataset,".permuted.tsv.gz",sep="")

df2<-fread(file2,header=TRUE,sep="\t")
df2<-df2[!is.na(df2$molecular_trait_id),]

# try first run without prunning
# # keep only nominal in permuted, note not all variants in permuted file available in nominal file
# vec<-df2$variant
# vec<-vec[!duplicated(vec)]
# 
# # note - if permuted file is not available do not truncate the nominal file afterwards
# to_prune<-df1$variant[which(df1$pvalue>0.5)]
# to_prune<-to_prune[which(!to_prune %in% vec)]
# 
# df1<-df1[which(!df1$variant %in% to_prune),]
# print(paste("N variants in permuted not in nominal:",nrow(df2[which(!df2$variant %in% df1$variant),]),sep=" "))
# 
# # exclude those as well - keep only those found in NON ge or microarray studies where nominal file is by default in V6 truncated
# df2<-df2[which(df2$variant %in% df1$variant),]

df1<-df1[,..columpsToKeep]

# make new directory if needed:
dir1<-paste(path,"nominal/",conditions,"/",sep="")
if (!dir.exists(dir1)) {
  system(paste("mkdir -p ",dir1,sep=""))
}

# save nominal file, and index:
file_out<-paste(path,"nominal/",conditions,"/",conditions,".tsv",sep="")
fwrite(df1,file=file_out,col.names = T,row.names = FALSE,sep="\t",quote = FALSE,na="NA")

# remove existing files
system(paste("rm ",file_out,".gz",sep=""))
system(paste("rm ",file_out,".gz.tbi",sep=""))

system(paste("bgzip ",file_out,sep=""))
system(paste("tabix -f -s13 -b14 -e14 -S1 ",file_out,".gz",sep=""))


#############################

df1$X<-paste(df1$molecular_trait_id,df1$variant,sep="_")
df2$X<-paste(df2$molecular_trait_id,df2$variant,sep="_")

df2<-df2[!duplicated(df2),]

print(paste("Number of tests in permutation not in original:",nrow(df2[which(!df2$X %in% df1$X),])))
print("N should be 0 if original file is not truncated")

# the qvalue is estimated in the whole gene set, so permuted file cannot be truncated, explore in the pipeline downstream if that is an issue
# df1<-df1[which(df1$X %in% df2$X),]


# should have been sorted above
# # deal with duplicated ids:
# dup<-df1$X[which(duplicated(df1$X))]
# dup<-dup[!duplicated(dup)]
# 
# df1$remove<-NA
# 
# # keep only first:
# for (i in 1:length(dup)) {
#   df1$remove[which(df1$X %in% dup[i])][2:length(df1$remove[which(df1$X %in% dup[i])])]<-1
# }
# 
# df1<-df1[which(is.na(df1$remove)),]

df<-merge(df1,df2,by="X",all.y=T)

print(paste("N variants in permuted: ",nrow(df2),sep=""))
print(paste("N variants in nominal available in permuted: ",nrow(df[!is.na(df$chromosome.x),]),sep=""))


rm(df1,df2)

# add TSS for each gene - same annotation as in eQTL catalog
gtf<-rtracklayer::import("/path/to/project")
gene<-as.data.frame(gtf)
rm(gtf)

gene$gene_id.2<-gsub("\\.[0-9]{1,2}$","",gene$gene_id)
gene<-gene[which( (gene$gene_id.2 %in% df$molecular_trait_object_id.y) & (gene$type=="gene")),c("gene_id.2","start","end","strand")]

gene$TSS<-NA
gene$TSS[which(gene$strand=="+")]<-gene$start[which(gene$strand=="+")]
gene$TSS[which(gene$strand=="-")]<-gene$end[which(gene$strand=="-")]

df<-merge(df,gene[,c("gene_id.2","TSS","strand")],by.x="molecular_trait_object_id.y",by.y="gene_id.2",all.x=T)

df$Distance_tss_variant<-df$position.y-df$TSS

df<-df[which(!is.na(df$strand)),]

# combine df and permuted data to replicate Nikos file, see how this file is used downstream by coloc pipeline
df_final<-data.frame(Phenotype_ID=df$molecular_trait_object_id.y,
                     Chromosome_phe=df$chromosome.y,
                     TSS=df$TSS,
                     TSS_end=df$TSS,
                     Strand=df$strand,
                     Total_no_variants_cis=df$n_variants,
                     Distance_tss_variant=df$Distance_tss_variant,
                     Best_variant_in_cis=df$variant.y,
                     Chromosome_var=df$chromosome.y,
                     Pos_variant=df$position.y,
                     Pos_variant_end=df$position.y,
                     DF=NA,
                     Dummy=df$variant.y,
                     Beta_dist_1=NA,
                     Beta_dist_2_number_of_ind_tests=NA,
                     Nominal_pvalue=df$pvalue.y,
                     Beta_regression=df$beta.y,
                     Empirical_pvalue=df$p_perm,
                     Corrected_pvalue=df$p_beta,
                     std.err=NA)


df_final<-df_final[which(!is.na(df_final$Corrected_pvalue) & df_final$Corrected_pvalue!=""),]

# # some corrected p-values set to 0, edit
# min_Corrected_pvalue<-min(min(df_final$Corrected_pvalue[which(df_final$Corrected_pvalue!=0)]),1E-322)
# df_final$Corrected_pvalue[which(df_final$Corrected_pvalue==0)]<-min_Corrected_pvalue

df_final$Best_variant_in_cis[which(df_final$Best_variant_in_cis==".")] <- as.character(df_final$Dummy[which(df_final$Best_variant_in_cis==".")])
df_final<-df_final[order(df_final$Chromosome_phe,df_final$Pos_variant,decreasing =F),]

print(paste("N FINAL variants in nominal available in permuted: ",nrow(df_final),sep=""))

# make new directory if needed:
dir2<-paste(path,"permuted/",sep="")
if (!dir.exists(dir2)) {
  system(paste("mkdir -p ",dir2,sep=""))
}

# save file, and index:
file_out_2<-paste(path,"permuted/",conditions,".permuted.txt",sep="")

fwrite(df_final,file_out_2,col.names = TRUE,row.names = FALSE,sep="\t",quote = FALSE,na="NA")

system(paste("rm ",file_out_2,".gz",sep=""))
system(paste("bgzip ",file_out_2))
rm(df_final)

####################

# for truncated files, save a key with links:

if (nrow(df[!is.na(df$chromosome.x),])!=nrow(df)) {
  
  df<-df[!is.na(df$chromosome.x),]
  
  df_final<-data.frame(Phenotype_ID=df$molecular_trait_object_id.y,
                       Chromosome_phe=df$chromosome.y,
                       TSS=df$TSS,
                       TSS_end=df$TSS,
                       Strand=df$strand,
                       Total_no_variants_cis=df$n_variants,
                       Distance_tss_variant=df$Distance_tss_variant,
                       Best_variant_in_cis=df$variant.y,
                       Chromosome_var=df$chromosome.y,
                       Pos_variant=df$position.y,
                       Pos_variant_end=df$position.y,
                       DF=NA,
                       Dummy=df$variant.y,
                       Beta_dist_1=NA,
                       Beta_dist_2_number_of_ind_tests=NA,
                       Nominal_pvalue=df$pvalue.y,
                       Beta_regression=df$beta.y,
                       Empirical_pvalue=df$p_perm,
                       Corrected_pvalue=df$p_beta,
                       std.err=NA)
  
  
  file_out_3<-paste(path,"permuted/",conditions,".permuted_in_nominal.txt",sep="")
  
  fwrite(df_final,file_out_3,col.names = TRUE,row.names = FALSE,sep="\t",quote = FALSE,na="NA")
  
  system(paste("rm ",file_out_3,".gz",sep=""))
  system(paste("bgzip ",file_out_3))
  
}

q("no")


