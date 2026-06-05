# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# MEM=35000
# bsub -Is -M"$MEM" -R"select[model==Intel_Platinum && mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group R \

library(data.table)
library(R.utils)
library(tidyr)

# only multithreaded when writting
setDTthreads(4)

rm(list=ls())

dataset<-"hu_2021"

path<-paste0("/path/to/project",dataset,"/")

df<-fread(paste(path,"raw_data/Merged.FDR.txt.gz",sep=""))

# create bed file to liftover - current poistion in b37/hg19
# A1 is the minor allele (effect allele) and the AllelFre is the minor allele frequency

df$variant<-paste("chr",df$Chr,"_",df$Pos,"_",df$Allele0,"_",df$Allele1,sep="")

bed<-df[,c("Chr","Pos","variant","Allele1","Allele0","AllelFre")]
bed<-bed[!duplicated(bed),]
bed$end<-bed$Pos+1

bed$Pos<-format(bed$Pos, scientific=F)
bed$end<-format(bed$end, scientific=F)

bed$chr<-paste("chr",bed$Chr,sep="")

setorder(bed,cols="Chr","Pos")

file_out<-paste(path,"other/",dataset,"_b37.bed",sep="")

write.table(bed[,c("chr","Pos","end","variant")],file_out,col.names=F,row.names=F,quote=F)


#### lift positions
system(paste("/software/team152/lf9/./liftOver ",
             file_out," /path/to/project ",
             file_out,"_lifted_hg38 ",file_out,"_no_lifted_hg38",sep=""))

system(paste("cut -f 4 ",file_out,"_no_lifted_hg38 | sed '/^#/d' > ",file_out,"_no_lifted_hg38_variants_to_exclude_tmp",sep=""))
system(paste("grep alt ",file_out,"_lifted_hg38 | cut -f 4 | cat - ",file_out,"_no_lifted_hg38_variants_to_exclude_tmp > ",
             file_out,"_no_lifted_hg38_variants_to_exclude",sep=""))

system(paste("wc -l ",file_out,"_no_lifted_hg38_variants_to_exclude",sep=""))
# 1559 ",path,"hu_2021//other/hu_2021_b37.bed_no_lifted_hg38_variants_to_exclude


### exclude non lifted:

bed_up<-fread(paste(file_out,"_lifted_hg38",sep=""),head=F)
colnames(bed_up)[2]<-"pos_b38"
colnames(bed_up)[4]<-"variant"

dim(bed)
# [1] 4188705       7
bed<-merge(bed,bed_up[,c("variant","pos_b38")],by="variant",all.y=T,sort=F)
dim(bed)
# [1] 4187769       8

# create tfam and tped like file to convert to vcf file
bed$morgan<-0
bed$genotype<-paste(bed$Allele1,bed$Allele0)

tped<-bed[,c("chr","variant","morgan","pos_b38","genotype")]

file_out_tped<-paste(path,"other/",dataset,"_b38.tped",sep="")
write.table(tped,file_out_tped,
            col.names=F,row.names=F,sep="\t",quote=F)
rm(tped)

tfam<-c("ID_1","ID_1","0","0","1","1")
tfam<-t(as.data.frame(tfam))
file_out_tfam<-paste(path,"other/",dataset,"_b38.tfam",sep="")
write.table(tfam,file_out_tfam,
            col.names=F,row.names=F,sep="\t",quote=F)
rm(tfam)

# create A1 allele
file_out_allele<-paste(path,"other/",dataset,"_b38_A1",sep="")
write.table(bed[,c("variant","Allele1")],file_out_allele,
            col.names=F,row.names=F,sep="\t",quote=F)


# create vcf:
system(paste("/software/team152/lf9/plink_linux_x86_64_20181202/./plink ",
             "--tfile ",path,"other/",dataset,"_b38 ",
             "--allow-no-sex ",
             "--a2-allele ",path,"other/",dataset,"_b38_A1 ",
             "--keep-allele-order --output-chr M --recode vcf-iid ",
             "--out ",path,"other/",dataset,"_b38",sep=""))


# # Note that most PLINK analyses treat the A1 (usually minor) allele as the reference allele, which makes sense when only biallelic variants are involved.
# # However, since it is conventional for VCF files to set the major allele as the reference allele instead
# 
# # # double check alleles, variants:
# export BCFTOOLS_PLUGINS=/software/team152/lf9/bcftools-1.16/plugins

system(paste0("bcftools +fixref ",path,"other/",dataset,"_b38.vcf ",
             "-Oz -o ",path,"other/",dataset,"_b38.vcf.gz -- -f /path/to/project -m top"))


# /software/team152/bcftools-1.16/./bcftools +fixref \
# ",path,"hu_2021/other/hu_2021_b38.vcf \
# -Oz -o ",path,"hu_2021/other/hu_2021_b38.vcf.gz \
# -- -f /path/to/project -m top
# # 
# # # SC, guessed strand convention
# # SC	TOP-compatible	0
# # SC	BOT-compatible	0
# # # ST, substitution types
# # ST	A>C	212611	5.1%
# # ST	A>G	935432	22.3%
# # ST	A>T	0	0.0%
# # ST	C>A	186858	4.5%
# # ST	C>G	0	0.0%
# # ST	C>T	760348	18.2%
# # ST	G>A	760856	18.2%
# # ST	G>C	0	0.0%
# # ST	G>T	186789	4.5%
# # ST	T>A	0	0.0%
# # ST	T>C	932208	22.3%
# # ST	T>G	212628	5.1%
# # # NS, Number of sites:
# # NS	total        	4187769
# # NS	ref match    	1054521	25.2%
# # NS	ref mismatch 	3133209	74.8%
# # NS	flipped      	1096	0.0%
# # NS	swapped      	3131045	74.8%
# # NS	flip+swap    	1068	0.0%
# # NS	unresolved   	0	0.0%
# # NS	fixed pos    	0	0.0%
# # NS	skipped      	39
# # NS	non-ACGT     	39
# # NS	non-SNP      	0
# # NS	non-biallelic	0
# 


# # Double check:
system(paste0("bcftools +fixref ",path,"other/",dataset,
"_b38.vcf.gz -- -f /path/to/project"))

# /software/team152/bcftools-1.16/./bcftools +fixref \
# ",path,"hu_2021/other/hu_2021_b38_2.vcf.gz \
# -- -f /path/to/project
# # # SC, guessed strand convention
# # SC	TOP-compatible	0
# # SC	BOT-compatible	0
# # # ST, substitution types
# # ST	A>C	189532	4.5%
# # ST	A>G	778624	18.6%
# # ST	A>T	0	0.0%
# # ST	C>A	209940	5.0%
# # ST	C>G	0	0.0%
# # ST	C>T	914976	21.8%
# # ST	G>A	917683	21.9%
# # ST	G>C	0	0.0%
# # ST	G>T	210031	5.0%
# # ST	T>A	0	0.0%
# # ST	T>C	777561	18.6%
# # ST	T>G	189383	4.5%
# # # NS, Number of sites:
# # NS	total        	4187730
# # NS	ref match    	4187730	100.0%
# # NS	ref mismatch 	0	0.0%
# # NS	skipped      	0
# # NS	non-ACGT     	0
# # NS	non-SNP      	0
# # NS	non-biallelic	0
# 
# 
# # #### VCF to BED
system(paste0("/software/team152/lf9/plink_linux_x86_64_20181202/./plink --vcf ",path,"other/",dataset,
"_b38.vcf.gz --keep-allele-order --allow-no-sex --double-id --make-bed --out ",path,"other/",dataset,"_b38_2"))

# /software/team152/plink_linux_x86_64_20181202/./plink \
# --vcf ",path,"hu_2021/other/hu_2021_b38_2.vcf.gz \
# --keep-allele-order --allow-no-sex \
# --double-id \
# --make-bed --out ",path,"hu_2021/other/hu_2021_b38_2


# # ##############################################################
# # # 16.2 UPDATE NAME VARIANTS TO CHR:POSITION_REF_ALT in b38

system(paste0("zcat ",path,"other/",dataset,"_b38.vcf.gz | cut -f '1-5' | awk '{print $3,$1\":\"$2\"_\"$4\"_\"$5}' > ",path,"other/",dataset,"_list_variants_b38"))

# # zcat ",path,"hu_2021/other/hu_2021_b38_2.vcf.gz | cut -f '1-5' | awk '{print $3,$1":"$2"_"$4"_"$5}' \
# # > ",path,"hu_2021/other/list_variants_hu_2021_b38

########################

bim<-fread(paste0(path,"other/",dataset,"_b38_2.bim"),head=F)
colnames(bim)[2]<-"variant"

ids<-fread(paste0(path,"other/",dataset,"_list_variants_b38"),head=F)
ids<-ids[-(1:2),]
colnames(ids)<-c("variant","variant_b38")

bim<-merge(bim,ids,by="variant",all.y=T)
bim<-bim[,c("variant","variant_b38")]

bed<-merge(bed,bim,by="variant",all.y=T)
rm(ids,bim)

### compare freq with freq in gnomad:

# system("gunzip /path/to/project")
system(paste0("awk 'NR==FNR{vals[$2];next} ($1) in vals' ",path,"other/",dataset,"_list_variants_b38 /path/to/project > ",path,"other/",dataset,"_list_variants_b38_gnomad"))
system(paste0("gzip ",path,"other/",dataset,"_list_variants_b38_gnomad"))

bed$Ref<-gsub("[0-9]{1,2}:[0-9]*_","",bed$variant_b38)
bed$Alt<-gsub("^[A-Z]*_","",bed$Ref)
bed$Ref<-gsub("_[A-Z]*$","",bed$Ref)

# no change - same freq
bed$Alt_freq[which(bed$Alt==bed$Allele1)]<-bed$AllelFre[which(bed$Alt==bed$Allele1)]

# flip - 1-freq
bed$Alt_freq[which(bed$Alt==bed$Allele0)]<-1-bed$AllelFre[which(bed$Alt==bed$Allele0)]

dim(bed[which(is.na(bed$Alt_freq)),])
# [1] 6382   14

# VARIANTS WHERE REF AND ALT HAVE BEEN SWAPPED -  same freq:

bed$Alt_freq[which( (bed$Allele0=="A") & (bed$Ref=="T") & (bed$Allele1=="C") & (bed$Alt=="G"))]<-bed$AllelFre[which( (bed$Allele0=="A") & (bed$Ref=="T") & (bed$Allele1=="C") & (bed$Alt=="G"))]
bed$Alt_freq[which( (bed$Allele0=="A") & (bed$Ref=="T") & (bed$Allele1=="G") & (bed$Alt=="C"))]<-bed$AllelFre[which( (bed$Allele0=="A") & (bed$Ref=="T") & (bed$Allele1=="G") & (bed$Alt=="C"))]
bed$Alt_freq[which( (bed$Allele0=="C") & (bed$Ref=="G") & (bed$Allele1=="A") & (bed$Alt=="T"))]<-bed$AllelFre[which( (bed$Allele0=="C") & (bed$Ref=="G") & (bed$Allele1=="A") & (bed$Alt=="T"))]
bed$Alt_freq[which( (bed$Allele0=="C") & (bed$Ref=="G") & (bed$Allele1=="T") & (bed$Alt=="A"))]<-bed$AllelFre[which( (bed$Allele0=="C") & (bed$Ref=="G") & (bed$Allele1=="T") & (bed$Alt=="A"))]
bed$Alt_freq[which( (bed$Allele0=="G") & (bed$Ref=="C") & (bed$Allele1=="T") & (bed$Alt=="A"))]<-bed$AllelFre[which( (bed$Allele0=="G") & (bed$Ref=="C") & (bed$Allele1=="T") & (bed$Alt=="A"))]
bed$Alt_freq[which( (bed$Allele0=="G") & (bed$Ref=="C") & (bed$Allele1=="A") & (bed$Alt=="T"))]<-bed$AllelFre[which( (bed$Allele0=="G") & (bed$Ref=="C") & (bed$Allele1=="A") & (bed$Alt=="T"))]
bed$Alt_freq[which( (bed$Allele0=="T") & (bed$Ref=="A") & (bed$Allele1=="C") & (bed$Alt=="G"))]<-bed$AllelFre[which( (bed$Allele0=="T") & (bed$Ref=="A") & (bed$Allele1=="C") & (bed$Alt=="G"))]
bed$Alt_freq[which( (bed$Allele0=="T") & (bed$Ref=="A") & (bed$Allele1=="G") & (bed$Alt=="C"))]<-bed$AllelFre[which( (bed$Allele0=="T") & (bed$Ref=="A") & (bed$Allele1=="G") & (bed$Alt=="C"))]

dim(bed[which(is.na(bed$Alt_freq)),])
# [1] 1096   14


# VARIANTS WHERE REF AND ALT HAVE BEEN SWAPPED AND FLIPPED - 1-freq:

bed$Alt_freq[which( (bed$Allele0=="A") & (bed$Ref=="G") & (bed$Allele1=="C") & (bed$Alt=="T"))]<-1-bed$AllelFre[which( (bed$Allele0=="A") & (bed$Ref=="G") & (bed$Allele1=="C") & (bed$Alt=="T"))]
bed$Alt_freq[which( (bed$Allele0=="A") & (bed$Ref=="C") & (bed$Allele1=="G") & (bed$Alt=="T"))]<-1-bed$AllelFre[which( (bed$Allele0=="A") & (bed$Ref=="C") & (bed$Allele1=="G") & (bed$Alt=="T"))]
bed$Alt_freq[which( (bed$Allele0=="C") & (bed$Ref=="T") & (bed$Allele1=="A") & (bed$Alt=="G"))]<-1-bed$AllelFre[which( (bed$Allele0=="C") & (bed$Ref=="T") & (bed$Allele1=="A") & (bed$Alt=="G"))]
bed$Alt_freq[which( (bed$Allele0=="C") & (bed$Ref=="A") & (bed$Allele1=="T") & (bed$Alt=="G"))]<-1-bed$AllelFre[which( (bed$Allele0=="C") & (bed$Ref=="A") & (bed$Allele1=="T") & (bed$Alt=="G"))]
bed$Alt_freq[which( (bed$Allele0=="G") & (bed$Ref=="A") & (bed$Allele1=="T") & (bed$Alt=="C"))]<-1-bed$AllelFre[which( (bed$Allele0=="G") & (bed$Ref=="A") & (bed$Allele1=="T") & (bed$Alt=="C"))]
bed$Alt_freq[which( (bed$Allele0=="G") & (bed$Ref=="T") & (bed$Allele1=="A") & (bed$Alt=="C"))]<-1-bed$AllelFre[which( (bed$Allele0=="G") & (bed$Ref=="T") & (bed$Allele1=="A") & (bed$Alt=="C"))]
bed$Alt_freq[which( (bed$Allele0=="T") & (bed$Ref=="G") & (bed$Allele1=="C") & (bed$Alt=="A"))]<-1-bed$AllelFre[which( (bed$Allele0=="T") & (bed$Ref=="G") & (bed$Allele1=="C") & (bed$Alt=="A"))]
bed$Alt_freq[which( (bed$Allele0=="T") & (bed$Ref=="C") & (bed$Allele1=="G") & (bed$Alt=="A"))]<-1-bed$AllelFre[which( (bed$Allele0=="T") & (bed$Ref=="C") & (bed$Allele1=="G") & (bed$Alt=="A"))]

dim(bed[which(is.na(bed$Alt_freq)),])
# [1]  0 15


gnomad<-fread(paste0(path,"other/",dataset,"_list_variants_b38_gnomad.gz"),head=F,sep=" ")
gnomad<-separate_wider_delim(gnomad, cols = "V1", delim = "\t", names = c("SNP", "CHROM"))
gnomad<-as.data.frame(gnomad)

colnames(gnomad)<-c("SNP","CHROM","POS","REF","ALT","AF","AF_nfe","AF_afr","AF_amr","AF_eas","AF_sas","AF_asj")
gnomad$AF_nfe<-as.numeric(as.character(gnomad$AF_nfe))

bed<-merge(bed,gnomad[,c("SNP","AF_nfe")],by.x="variant_b38",by.y="SNP",all.x=T,sort=F)

bed$value<-((bed$Alt_freq-bed$AF_nfe)^2)/((bed$Alt_freq+bed$AF_nfe)*(2-bed$Alt_freq-bed$AF_nfe))

summary(bed$value)

cbPalette <- c("#999999", "#E69F00", "#56B4E9")
bed$col[which(bed$value<0.125)]<-0
bed$col[which(bed$value>=0.125)]<-1

table(bed$col,useNA="ifany")
# 

bed$col<-as.character(bed$col)

library(ggplot2)
p3<-ggplot(bed, aes(y=Alt_freq, x=AF_nfe)) +
  geom_point(aes(colour = col)) + ylab(paste0(dataset," FRQ Alt")) + xlab("Gnomad FRQ Alt") + scale_colour_manual(values=cbPalette)

pdf(paste0(path,"other/plot_alt_freq_",dataset,"_Gnomad.pdf"),width = 6, height = 5)
print(p3)
dev.off()

fwrite(bed,paste0(path,"other/snp_info_b37_b38_tmp.txt.gz"),
       col.names=T,row.names=F,sep="\t",quote=F)

q("no")
