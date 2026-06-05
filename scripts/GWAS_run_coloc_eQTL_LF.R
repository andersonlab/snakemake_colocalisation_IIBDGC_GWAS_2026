# Adapted by: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# debugging - modules to load:
# 
# # singularity exec iibdgc_postprocess_10_singularity.sif
# # singularity exec iibdgc_postprocess_10_singularity.sif
# # singularity exec iibdgc_postprocess_10_singularity.sif

# # debugging - launch interacive job:
# MEM=15000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R


suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("purrr"))
suppressPackageStartupMessages(library("coloc"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("devtools"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("purrrlyr"))
suppressPackageStartupMessages(library("Rsamtools"))



#Parse command-line options
option_list <- list(
  make_option(c("-p", "--phenotype"), type="character", default=NULL,
              help="Type of QTLs used for coloc.", metavar = "type"),
  make_option(c("-w", "--window"), type="character", default=NULL,
              help="Size of the cis window.", metavar = "type"),
  make_option(c("--gwas"), type="character", default=NULL,
              help="Name of the GWAS trait", metavar = "type"),
  make_option(c("-d", "--dir"), type="character", default=NULL,
              help="Path to GWAS summary stats directory.", metavar = "type"),
  make_option(c("--qtl"), type="character", default=NULL,
              help="Path to the QTL directory.", metavar = "type"),
  make_option(c("-o", "--outdir"), type="character", default=NULL,
              help="Path to the output directory.", metavar = "type"),
  make_option(c("-s", "--samplesizes"), type="character", default=NULL,
              help="Path to the tab-separated text file with condition names and sample sizes.", metavar = "type"),
  make_option(c( "--gwasvarinfo"), type="character", default=NULL,
              help="Variant infromation file for the GWAS dataset.", metavar = "type"),
  make_option(c("--qtlvarinfo"), type="character", default=NULL,
              help="Variant information file for the QTL dataset.", metavar = "type"),
  make_option(c("--gwaslist"), type="character", default=NULL,
              help="Path to the list of GWAS studies.", metavar = "type"),
  make_option(c("--chunk"), type="integer", default=NULL,
              help="chunk to run", metavar = "type"),
  make_option(c("--numberGWAS"), type="integer", default=NULL,
              help="numberGWAS", metavar = "type"),
  make_option(c("--function_path_source"), type="character", default=NULL,
              help="Path to the functions_me_GTEx file", metavar = "type")

)
opt <- parse_args(OptionParser(option_list=option_list))

#Debugging - example sc-eQTL data with <20 chunks 
#opt = list(gwas = "CD", w = "2e6", p = "featureCounts", d = "/path/to/project",
          #  o = "/path/to/project",
          #  qtl = "/path/to/project",
          #  s = "/path/to/project",#
          #  gwasvarinfo = "/path/to/project",
          #  qtlvarinfo = "/path/to/project",
          #  gwaslist = "/path/to/project",
          #  chunk=1,function_path_source="/path/to/project")


# opt = list(gwas = "CD", w = "2e6", p = "featureCounts", 
#            d = "/path/to/project",
#            o = "/path/to/project",
#            qtl = "/path/to/project",
#            s = "/path/to/project",#
#            gwasvarinfo = "/path/to/project",
#            qtlvarinfo = "/path/to/project",
#            gwaslist = "/path/to/project",
#            chunk=20,function_path_source="/path/to/project")

# opt = list(gwas = "CD", w = "2e6", p = "featureCounts", 
#            d = "/path/to/project",
#            o = "/path/to/project",
#            qtl = "/path/to/project",
#            s = "/path/to/project",#
#            gwasvarinfo = "/path/to/project",
#            qtlvarinfo = "/path/to/project",
#            gwaslist = "/path/to/project",
#            chunk=7,function_path_source="/path/to/project")

# opt = list(gwas = "CD", w = "2e6", p = "featureCounts", 
#            d = "/path/to/project",
#            o = "/path/to/project",
#            qtl = "/path/to/project",
#            s = "/path/to/project",#
#            gwasvarinfo = "/path/to/project",
#            qtlvarinfo = "/path/to/project",
#            gwaslist = "/path/to/project",
#            chunk=19,function_path_source="/path/to/project")

 #opt$d = "/path/to/project"


#Extract parameters for CMD options
gwas_id = opt$gwas
cis_window = as.numeric(opt$w)
phenotype = opt$p
gwas_dir = opt$d
qtl_dir = opt$qtl
outdir = opt$o
sample_size_path = opt$s
gwas_var_path = opt$gwasvarinfo
qtl_var_path = opt$qtlvarinfo
gwas_list = opt$gwaslist
chunk=opt$chunk
numberGWAS=opt$numberGWAS
function_path=opt$function_path_source

source(function_path)

#Import variant information
#gwas_var_info = importVariantInformation(gwas_var_path) # This set of variants is with v37 of the genome
qtl_var_info = importVariantInformation(qtl_var_path) # v38
print(head(qtl_var_info))
print("Variant information imported.")
#macromap_eQTLs<-read.table("~/myscratch/MacroMap/Analysis/eQTLs/Macromap_fds/analysis/mashR/files/8_flashR_condition_by_condition/rownames_sign_eQTLs_betas.txt",h=F)
#macromap_eQTLs_genes<- macromap_eQTLs_genes<-sapply(strsplit(unique(substr(macromap_eQTLs$V1,1,15)),"[,]"),"[[",1)

#Import list of GWAS studies
gwas_stats_labeled = readr::read_tsv(gwas_list, col_names = c("trait","file_name","type"), col_type = "ccc")

#Import sample sizes
sample_sizes = readr::read_tsv(sample_size_path, col_names = c("condition_name", "sample_size"), col_types = "ci")
sample_sizes_list = as.list(sample_sizes$sample_size)
names(sample_sizes_list) = sample_sizes$condition_name
sample_sizes_list = sample_sizes_list[1]

#Construct a new QTL list
phenotype_values = constructQtlListForColoc(qtl_dir, sample_sizes_list)
print(head(phenotype_values))

#Spcecify the location of the GWAS summary stats file
gwas_file_name = dplyr::filter(gwas_stats_labeled, trait == gwas_id)$file_name
gwas_prefix = file.path(gwas_dir, gwas_file_name)

#Prefilter coloc candidates
qtl_df_list = prefilterColocCandidates(phenotype_values$min_pvalues, gwas_prefix,
                                       gwas_variant_info = qtl_var_info, fdr_thresh = 0.05,
                                       overlap_dist = 1e7, gwas_thresh = 1e-5)
qtl_pairs = purrr::map_df(qtl_df_list, identity) %>% unique()
names(qtl_pairs)<-c("phenotype_id","snp_id")
print(head(qtl_pairs))

# chunking
#qtl_pairs<-qtl_pairs[qtl_pairs$phenotype_id %in% macromap_eQTLs_genes,]
d<-1:length(qtl_pairs$phenotype_id)
d.split<-split(d, sort(d%%20))

# for sc-eqtl data, create an empty output file if now genes in the chunck:

if (length(d.split)<chunk) {

    print(paste("Pre-filtering completed, no genes to test in chunk",chunk))
    print("Coloc completed.")
    vec<-c("condition_name","phenotype_id","snp_id",".row","nsnps","PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf","qtl_pval","gwas_pval","qtl_lead","gwas_lead","chr","gwas_lead_pos","qtl_lead_pos","gwas_trait")
    coloc_hits<-as.data.frame(matrix(nrow=0,ncol=length(vec)))
    colnames(coloc_hits)<-vec
    coloc_output = file.path(outdir, paste(gwas_id, phenotype, opt$w,sample_sizes$condition_name,"chunk",chunk,"txt", sep = "."))
    write.table(coloc_hits, coloc_output, sep = "\t", quote = FALSE, row.names = FALSE)

} else {

    qtl_pairs<-qtl_pairs[d.split[[chunk]],]

    print(paste("Pre-filtering completed, number of genes to test",length(d.split[[chunk]]),"chunk",chunk))

    options(warnings=1) 
    #Test for coloc
    coloc_res_list = purrr::map2(phenotype_values$qtl_summary_list, phenotype_values$sample_sizes,
                                ~colocMolecularQTLsByRow_nikos(qtl_pairs, qtl_summary_path = .x,
                                                                gwas_summary_path = paste0(gwas_prefix, ".txt.gz"),
                                                                gwas_variant_info = qtl_var_info,
                                                                qtl_variant_info = qtl_var_info,
                                                                N_qtl = .y, cis_dist = cis_window,N=numberGWAS))

    print("Coloc completed.")

    #Export results
    coloc_hits = purrr::map_df(coloc_res_list, identity, .id = "condition_name") %>% dplyr::arrange(gwas_lead) %>% dplyr::mutate(gwas_trait=gwas_id)
    coloc_output = file.path(outdir, paste(gwas_id, phenotype, opt$w,sample_sizes$condition_name,"chunk",chunk,"txt", sep = "."))
    write.table(coloc_hits, coloc_output, sep = "\t", quote = FALSE, row.names = FALSE)

}
