'''
converts arms files from mips into psl alignments that can then be visualized on
a genome browser (e.g. igv or ucsc genome browser)

Reports genomic coordinates of individual mips for later use
'''

configfile: 'map_mip_arms.yaml'
out_dir=config['output_folder']

rule all:
	input:
		filtered_bed=out_dir+'/'+config['panel_name']+'_'+config['genome_name']+'_mips_filtered.bed',
		hub_file=out_dir+'/'+config['panel_name']+'_'+config['genome_name']+'_hub.txt',
		snakemake_file=out_dir+'/snakemake_params/map_mip_arms.smk'

rule copy_files:
	'''
	copies snakemake parameters to output folder for reproducibility
	'''
	input:
		yaml_file='map_mip_arms.yaml',
		snakemake_file='map_mip_arms.smk',
		scripts_folder='scripts',
		env_folder='envs'
	output:
		yaml_file=out_dir+'/snakemake_params/map_mip_arms.yaml',
		snakemake_file=out_dir+'/snakemake_params/map_mip_arms.smk',
		scripts_folder=directory(out_dir+'/snakemake_params/scripts'),
		env_folder=directory(out_dir+'/snakemake_params/envs')
	shell:
		'''
		cp {input.yaml_file} {output.yaml_file}
		cp {input.snakemake_file} {output.snakemake_file}
		cp -R {input.scripts_folder} {output.scripts_folder}
		cp -R {input.env_folder} {output.env_folder}
		'''


rule arms_to_fasta:
	'''
	converts ligation and extension arms into fasta sequences for mapping onto
	the genome. Reverse complements the ligation arm and concatenates the two
	arms as 'exons'
	'''
	input:
		arms_file=config['mip_arms_file']
	output:
		mips_fasta=out_dir+'/'+config['panel_name']+'_mips.fa'
	script:
		'scripts/arms_to_fasta.py'

rule map_arms_to_genome:
	'''
	maps the fasta version of mip arms to the genome using blat
	'''
	input:
		mips_fasta=out_dir+'/'+config['panel_name']+'_mips.fa',
		genome_path=config['genome_path']
	output:
		mips_psl=out_dir+'/'+config['panel_name']+'_'+config['genome_name']+'_mips.psl'
	shell:
		'blat {input.genome_path} {input.mips_fasta} -minScore=5 -tileSize=6 {output.mips_psl}'

rule filter_best_hit:
	'''
	searches the blat hits and gets the best blat hit (or best tied hits).
	Also creates a coordinates file with the genomic coordinates of every mip
	'''
	input:
		mips_psl=out_dir+'/'+config['panel_name']+'_'+config['genome_name']+'_mips.psl'
	params:
		error_count=config['error_count']
	output:
		filtered_psl=out_dir+'/'+config['panel_name']+'_'+config['genome_name']+'_mips_filtered.psl',
		coords_file=out_dir+'/'+config['panel_name']+'_'+config['genome_name']+'_genomic_coords.tsv'
	script:
		'scripts/filter_best_hit.py'

rule psl_to_bed:
	'''
	converts a psl file into a bed file for other downstream analysis programs
	'''
	input:
		filtered_psl=out_dir+'/'+config['panel_name']+'_'+config['genome_name']+'_mips_filtered.psl'
	conda:
		'envs/pslToBed.yaml'
	output:
		filtered_bed=out_dir+'/'+config['panel_name']+'_'+config['genome_name']+'_mips_filtered.bed'
	shell:
		'pslToBed {input.filtered_psl} {output.filtered_bed}'

rule psl_to_big_psl:
	input:
		filtered_psl=out_dir+'/'+config['panel_name']+'_'+config['genome_name']+'_mips_filtered.psl'
	conda:
		'envs/pslToBigPsl.yaml'
	output:
		big_psl=out_dir+'/'+config['panel_name']+'_'+config['genome_name']+'_mips_filtered.bigPsl.txt'
	shell:
		'pslToBigPsl {input.filtered_psl} stdout | sort -k1,1 -k2,2n > {output.big_psl}'

rule bed_to_big_bed:
	input:
		big_psl=out_dir+'/'+config['panel_name']+'_'+config['genome_name']+'_mips_filtered.bigPsl.txt'
	conda:
		'envs/bedToBigBed.yaml'
	output:
		big_bed=out_dir+'/'+config['panel_name']+'_'+config['genome_name']+'_mips_filtered.bigPsl.bb'
	shell:
		'bedToBigBed -as=bigPsl.as -type=bed12+13 -tab {input.big_psl} chrom.sizes {output.big_bed}'

rule make_hub:
	'''
	creates a hub.txt file for viewing data on the UCSC genome browser
	'''
	input:
		big_bed=out_dir+'/'+config['panel_name']+'_'+config['genome_name']+'_mips_filtered.bigPsl.bb'
	params:
		hub_name=config['panel_name'],
		email=config['email'],
		genome=config['genome_name'],
		big_bed=config['panel_name']+'_'+config['genome_name']+'_mips_filtered.bigPsl.bb'
	output:
		hub_file=out_dir+'/'+config['panel_name']+'_'+config['genome_name']+'_hub.txt'
	script:
		'scripts/make_hub.py'
