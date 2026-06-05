# get input dir from config
import yaml
import os

configfile_yaml = open("/path/to/project").read()
config = yaml.load(configfile_yaml, Loader=yaml.FullLoader)

eqtl_study = config['eqtl_study']
gwas_study = config['gwas_study']

# general projects
input_dir = os.path.join(config['input_dir'], eqtl_study)

output_dir = os.path.join(config["input_dir"],"results", eqtl_study+"_"+gwas_study,"coloc_results")
script_dir = os.path.join(config["input_dir"],"results", eqtl_study+"_"+gwas_study)
gwas_dir = os.path.join(config["input_dir"],config["gwas_dir"])


# managed access projects - uncomment the following line:
# input_dir = os.path.join(config['input_dir'], "QTL_managed_access", eqtl_study)

condition_file = os.path.join(input_dir, config["condition_file"])

print(condition_file)

# Get the condition list
with open(condition_file, "r") as f:
    conditions = [line.strip() for line in f if line.strip()]

print(len(conditions))

# wildcards as definitions where the value will iterate over, defined as overall rules:

localrules: run_all
rule run_all:
	input:
		expand("{output_dir}/{gwas}.{phenotype}.{coloc_window}.{condition}.chunk.{chunks}.txt", gwas = config["gwas_traits"], phenotype = config["coloc_phenotypes"],  coloc_window = config["coloc_window"], output_dir=output_dir,condition = conditions, chunks = config["chunks"] )
	output:
		output_dir + "/coloc.DONE"
	resources:
		mem = 100
	threads: 1
	shell:
		"echo 'Done!' > {output}"



rule run_coloc:
	input:
		input_dir + "/nominal"
	output:
		"{output_dir}/{gwas}.{phenotype}.{coloc_window}.{condition}.chunk.{chunks}.txt"
	params:
		output_dir=output_dir,
		input_dir=input_dir,
		gwas_dir=gwas_dir,
		phenotype=config["coloc_phenotypes"],
		coloc_window=config["coloc_window"],
		script_dir=script_dir
	resources:
		queue="normal",
		mem_mb=15000,
		mem=15000, 
		mem_mib=15000,
		disk_mb=15000,
		tmpdir="tmp"
	threads: 1
	singularity:
		"/software/hgi/softpack/installs/groups/team152/snakemake_coloc_no_sm/4-scripts/singularity.sif"
	shell:
		r"""
			Rscript {params.script_dir}/GWAS_run_coloc_eQTL_LF.R \
				--phenotype {params.phenotype} --window {wildcards.coloc_window} \
				--gwas {wildcards.gwas} --dir {params.gwas_dir} --outdir {params.output_dir} \
				--qtl {params.input_dir} --samplesizes {params.input_dir}/sample_size_per_condition/{wildcards.condition} \
				--gwasvarinfo {params.input_dir}/variant_info/{wildcards.condition}.variant.info.tsv.gz --qtlvarinfo {params.input_dir}/variant_info/{wildcards.condition}.variant.info.tsv.gz \
				--gwaslist {params.gwas_dir}/{config[gwas_list]} --chunk {wildcards.chunks} --function_path_source {params.script_dir}/functions_me_eQTL_bh.R
		"""