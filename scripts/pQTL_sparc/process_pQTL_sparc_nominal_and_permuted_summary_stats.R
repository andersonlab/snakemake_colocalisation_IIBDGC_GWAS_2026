# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#


# # singularity exec iibdgc_postprocess_10_singularity.sif

# MEM=25000
# bsub -Is -M"$MEM" -R"select[model==Intel_Platinum && mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q long R

library(data.table)
library(R.utils)
library(dplyr)
library(rtracklayer)
library(qvalue)
library(stringr)

rm(list=ls())


# provide trait (in this case decode or ukb pqtl data)
args = commandArgs(trailingOnly=TRUE)

dataset<-args[1]
print(dataset)

condition<-args[1]
print(condition)

# dataset<-"pQTL_sparc"
# condition<-"plasma_ibd_patients"


path_gwas<-"/path/to/project"
path<-"/path/to/project"
path_proteomics<-paste0(path,dataset,"/")

# use gencode to add TSS start site:
gtf<-rtracklayer::import(paste(path_gwas,'post_imputation/2022/analysis/metaanalysis/annotation/gencode.v39.annotation.gtf.gz',sep=""))
gene<-as.data.frame(gtf)
rm(gtf)

gene<-gene[which(gene$type=="gene"),]
gene<-gene[which(!duplicated(gene$gene_name)),]

gene$TSS<-gene$start
gene$TSS[which(gene$strand=="-")]<-gene$end[which(gene$strand=="-")]


# get list of files to load:
pt_list<-fread("/path/to/project",head=F)
dim(pt_list[which(pt_list$V1 %in% gene$gene_name),])
# 2894
dim(pt_list[which(!pt_list$V1 %in% gene$gene_name),])
# [1] 27  1

pt_list[which(!pt_list$V1 %in% gene$gene_name),]
#                    V1
#                <char>
#  1:         FUT3_FUT5
#  2:             GATD3
#  3:            GPR15L
#  4:       IL12A_IL12B
#  5:           KIR2DL2
#  6:              LEG1
#  7:    LGALS7_LGALS7B
#  8:            LILRA3
#  9:              MENT
# 10:         MICB_MICA
# 11:          NTproBNP
# 12:             PALM2
# 13:              SARG
# 14:    SPACA5_SPACA5B
# 15:              WARS
# 16: AMY1A_AMY1B_AMY1C
# 17:             BAP18
# 18:      BOLA2_BOLA2B
# 19:              CERT
# 20:    CGB3_CGB5_CGB8
# 21:     CKMT1A_CKMT1B
# 22:     CTAG1A_CTAG1B
# 23:      DEFA1_DEFA1B
# 24: DEFB103A_DEFB103B
# 25: DEFB104A_DEFB104B
# 26:     DEFB4A_DEFB4B
# 27:         EBI3_IL27
#                    V1

tmp<-pt_list[which(!pt_list$V1 %in% gene$gene_name),]
tmp<-tmp[grep("_",tmp$V1),]
tmp<-tmp[order(tmp$V1),]

tmp2<-unlist(strsplit(tmp$V1,"_"))

tmp<-rbind(tmp,tmp)
tmp<-rbind(tmp,as.data.frame(as.matrix(c("CGB3_CGB5_CGB8","AMY1A_AMY1B_AMY1C"))))


tmp<-tmp[order(tmp$V1),]
tmp$gene_name<-tmp2

pt_list<-pt_list[which(!pt_list$V1 %in% tmp$V1),]
pt_list$gene_name<-pt_list$V1

pt_list<-rbind(pt_list,tmp)

# GATD3 = GATD3A + GATD3B
# GPR15L = GPR15LG
# SARG = C1orf116
# BAP18 = BACC1

pt_list$gene_name[which(pt_list$V1=="GPR15L")]<-"GPR15LG"
pt_list$gene_name[which(pt_list$V1=="SARG")]<-"C1orf116"
pt_list$gene_name[which(pt_list$V1=="BAP18")]<-"BACC1"

pt_list$gene_name[which(pt_list$V1=="GATD3")]<-"GATD3A"
tmp<-as.data.frame(as.matrix(t(c("GATD3","GATD3B"))))
colnames(tmp)<-colnames(pt_list)
pt_list<-rbind(pt_list,tmp)

colnames(pt_list)[1]<-"protein_name"

pt_list<-pt_list[!duplicated(pt_list$gene_name),]


gene<-gene[which(gene$gene_name %in% pt_list$gene_name),]
gene<-gene[which(gene$type=="gene"),]

dim(gene)
gene<-merge(pt_list,gene,by="gene_name")
dim(gene)

gene$TSS<-gene$start
gene$TSS[which(gene$strand=="-")]<-gene$end[which(gene$strand=="-")]

table(gene$seqnames)
#  chr1  chr2  chr3  chr4  chr5  chr6  chr7  chr8  chr9 chr10 chr11 chr12 chr13 
#   337   196   137   127   128   155   129   110   119   111   179   153    48 
# chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22  chrX  chrY  chrM 
#    91    75   119   168    34   238    77    30    75    92     0     0 

# only retain autosomal data - pQTL only done for those:
gene<-gene[which(!gene$seqnames %in% c("chrX","chrY","chrM")),]


# retrieve cis data for all genes:
for (i in 1:nrow(gene)) {

    # get the list of files per protein
    list_files<-list.files(paste0(path_proteomics,"raw/pQTLstatsForLaura/",gene$protein_name[i],"/"))

    # retain only the chr where the gene is:
    list_files<-list_files[grep(paste0("_",gene$seqnames[i],"_"),list_files)]

    file_name<-paste0(path_proteomics,"raw/pQTLstatsForLaura/",gene$protein_name[i],"/",list_files)
    print(file_name)

    if (file.exists(file_name)) {

        tmp<-fread(file_name)
        tmp<-tmp[which(tmp$GENPOS>=gene$TSS[i]-1000000 & tmp$GENPOS<=gene$TSS[i]+1000000),]

        if (nrow(tmp)>0) {
          tmp$gene_id<-gene$gene_name[i]
          tmp$Pval<-10^-tmp$LOG10P
          
          tmp<-tmp[which(!is.na(tmp$Pval)),]
          tmp$Pval[which(tmp$Pval==0)]<-1E-320
          
          tmp$qval<-qvalue(tmp$Pval)$qval

          tmp$TSS<-gene$TSS[i]
          tmp$strand<-gene$strand[i]
        }
        
    }

    if (i==1) {
        d<-tmp
    } else {
        d<-rbind(d,tmp,fill=T)
    }
    rm(file_name,tmp,list_files)


}


d$variant_b38<-paste0(d$CHROM,"_",d$GENPOS,"_",d$ALLELE0,"_",d$ALLELE1)


d$an<-2*d$N
d$ac<-round((d$an)*(d$A1FREQ))

d$r2<-NA
d$pvalue<-d$Pval
d$molecular_trait_object_id<-d$gene_id
d$molecular_trait_id<-d$gene_id
d$maf<-pmin(d$A1FREQ,1-d$A1FREQ)

d$gene_id<-d$gene_id
d$median_tpm<-NA
d$beta<-d$BETA # beta for effect allele, which is the reference (although SNP name is chr_pos_Alt_Ref and chr_pos_effect_other)
d$se<-d$SE
d$chromosome<-d$CHROM
d$position<-d$GENPOS
d$ref<-d$ALLELE0
d$alt<-d$ALLELE1

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
d$chromosome<-as.numeric(gsub("chr","",d$CHROM))


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
# [1] 2830

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
#[1] 2830

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
#[1] "N FINAL variants in nominal available in permuted: 2830"

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
dir4<-paste(path,dataset,"/sample_size_per_condition/",sep="")
if (!dir.exists(dir4)) {
  system(paste("mkdir -p ",dir4,sep=""))
}

file.out<-paste(path,dataset,"/sample_size_per_condition/",d3$dataset,sep="")
fwrite(d3,file.out,col.names=F,row.names=F,quote=F,sep="\t")

####################

q("no")


