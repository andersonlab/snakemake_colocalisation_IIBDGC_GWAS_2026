# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# # singularity exec iibdgc_postprocess_10_singularity.sif

# MEM=15000
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

# use gencode to add TSS start site:
gtf<-rtracklayer::import(paste(path_gwas,'post_imputation/2022/analysis/metaanalysis/annotation/gencode.v39.annotation.gtf.gz',sep=""))
gene<-as.data.frame(gtf)
rm(gtf)

gene<-gene[which(gene$type=="gene"),]
gene<-gene[which(!duplicated(gene$gene_name)),]

gene$TSS<-gene$start
gene$TSS[which(gene$strand=="-")]<-gene$end[which(gene$strand=="-")]


# get list of files to load:
pt<-list.files(paste0(path_proteomics,condition,"/"))
pt<-as.data.frame(pt[grep(".txt.gz",pt)])
colnames(pt)<-c("files")
nrow(pt)

pt$gene_symbol<-gsub("GBR_UKB_OLINK2_OID[0-9]*_","",pt$files)
pt$gene_symbol<-gsub("GBR_UKB_OLINK_OID[0-9]*_","",pt$gene_symbol)
pt$gene_symbol<-gsub("Proteomics_SMP_PC0_[0-9]*_[0-9]*_","",pt$gene_symbol)
pt$gene_symbol<-gsub("_adjAgeSexBatPC_InvNorm_22122022.txt.gz","",pt$gene_symbol)
pt$gene_symbol<-gsub("_10032022.txt.gz","",pt$gene_symbol)
pt$gene_name1<-gsub("_.*","",pt$gene_symbol)
pt$gene_name2<-gsub(".*_","",pt$gene_symbol)


# pt<-pt[which(pt$gene_symbol!="NA"),]
# dim(pt[which(pt$gene_symbol=="NA"),])
# # 0

pt$source<-NA
pt$source[grep("GBR_UKB_",pt$files)]<-"ukb"
pt$source[grep("Proteomics_SMP_PC0_",pt$files)]<-dataset
table(pt$source,useNA="ifany")
#  ukb 
# 2931 
# decode 
#   5235 


dim(pt[which(!pt$gene_name1 %in% gene$gene_name),])
# [1] 15  3
# [1] 245   3
dim(pt[which(pt$gene_name1 %in% gene$gene_name),])
# [1] 2916    3
# [1] 4990    3

pt1<-merge(gene[,c("gene_name","gene_id","start","end","strand","seqnames","TSS","gene_id")],pt[,c("files","gene_name1")],by.x="gene_name",by.y="gene_name1")
pt2<-merge(gene[,c("gene_name","gene_id","start","end","strand","seqnames","TSS","gene_id")],pt[,c("files","gene_name2")],by.x="gene_name",by.y="gene_name2")

pt<-rbind(pt1,pt2)
pt<-pt[!duplicated(pt),]
dim(pt)
# [1] 2918    9

# exclude genes in chrY:
pt<-pt[which(!pt$seqnames %in% c("chrY","chrM","")),]

# i=1359
# i=1326

# rename the files to load bgz - created by HGI
pt$files<-gsub(".txt.gz",".txt.bgz",pt$files)

# retrieve cis data for all genes:
for (i in 1:nrow(pt)) {

  file_name<-paste0(path_proteomics,condition,"/",pt$files[i])

  if (file.exists(file_name)) {
    
    # tmp<-fread(file_name)
    # tmp<-tmp[which(tmp$Chrom==pt$seqnames[i] & tmp$Pos>=pt$TSS[i]-1000000 & tmp$Pos<=pt$TSS[i]+1000000),]


    param<-GRanges(c(pt$seqnames[i]), IRanges((pt$TSS[i]-1000000):(pt$TSS[i]+1000000)))
    tbx<-Rsamtools::TabixFile(file_name)
        
    res <- Rsamtools::scanTabix(tbx, param=param)

    tmp <- Map(function(elt) {
        read.csv(textConnection(elt), sep="\t", header=T)
        }, res)
    tmp<-as.data.table(tmp[1])
    colnames(tmp)<-c("Chrom","Pos","Name","rsids","effectAllele","otherAllele","Beta","Pval","minus_log10_pval","SE","N","ImpMAF")


    if (nrow(tmp)>0) {
      tmp$gene_id<-as.character(pt$gene_name[i])
      tmp<-tmp[which(!is.na(tmp$Pval)),]
      tmp$Pval[which(tmp$Pval==0)]<-1E-320
      tmp<-tmp[which(!is.na(tmp$Pval)),]
      tmp$qval<-qvalue(tmp$Pval)$qval
      tmp$TSS<-pt$TSS[i]
      tmp$strand<-as.character(pt$strand[i])
    }

    if (i==1) {
      d<-tmp
    } else {
      d<-rbind(d,tmp)
    }
    rm(tmp)
  } else {
    print("File",file_name,"not found")
  }

  rm(file_name)

}

d$variant_b38<-paste0(d$Chrom,"_",d$Pos,"_",d$otherAllele,"_",d$effectAllele)


d$an<-2*d$N
d$ac<-round((d$an)*(d$ImpMAF))

d$r2<-NA
d$pvalue<-d$Pval
d$molecular_trait_object_id<-d$gene_id
d$molecular_trait_id<-d$gene_id
d$maf<-d$ImpMAF
d$gene_id<-d$gene_id
d$median_tpm<-NA
d$beta<-d$Beta # beta for effect allele, which is the reference (although SNP name is chr_pos_Alt_Ref and chr_pos_effect_other)
d$se<-d$SE
d$chromosome<-d$Chrom
d$position<-d$Pos
d$ref<-d$otherAllele
d$alt<-d$effectAllele

d$type<-"indel"
d$type[which(nchar(d$ref)==1 & nchar(d$alt)==1)]<-"SNP"

d$rsid<-d$variant_b38
d$variant<-d$variant_b38

rm(gene)

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
dir4<-paste(path,dataset,"/sample_size_per_tissue/",sep="")
if (!dir.exists(dir4)) {
  system(paste("mkdir -p ",dir4,sep=""))
}

file.out<-paste(path,dataset,"/sample_size_per_tissue/",d3$dataset,sep="")
fwrite(d3,file.out,col.names=F,row.names=F,quote=F,sep="\t")

####################

q("no")


