'''
converts arms files from mips into psl alignments that can then be visualized on
a genome browser (e.g. igv or ucsc genome browser)

Reports genomic coordinates of individual mips for later use
'''

configfile: 'map_mip_arms.yaml'
output=config['output_folder']

rule all:
	input:
		filtered_psl=output+'/'+config['panel_name']+'_'+config['genome_name']+'_mips_filtered.psl'
rule arms_to_fasta:
	'''
	converts ligation and extension arms into fasta sequences for mapping onto
	the genome. Reverse complements the ligation arm and concatenates the two
	arms as 'exons'
	'''
	input:
		arms_file=config['mip_arms_file']
	output:
		mips_fasta=output+'/'+config['panel_name']+'_mips.fa'
	script:
		'scripts/arms_to_fasta.py'

rule map_arms_to_genome:
	'''
	maps the fasta version of mip arms to the genome using blat
	'''
	input:
		mips_fasta=output+'/'+config['panel_name']+'_mips.fa',
		genome_path=config['genome_path']
	output:
		mips_psl=output+'/'+config['panel_name']+'_'+config['genome_name']+'_mips.psl'
	shell:
		'blat {input.genome_path} {input.mips_fasta} -minScore=5 -tileSize=6 {output.mips_psl}'

rule filter_best_hit:
	'''
	searches the blat hits and gets the best blat hit (or best tied hits).
	Also creates a coordinates file with the genomic coordinates of every mip
	'''
	input:
		mips_psl=output+'/'+config['panel_name']+'_'+config['genome_name']+'_mips.psl'
	output:
		filtered_psl=output+'/'+config['panel_name']+'_'+config['genome_name']+'_mips_filtered.psl',
		coords_file=output+'/'+config['panel_name']+'_'+config['genome_name']+'_genomic_coords.tsv'
	script:
		'scripts/filter_best_hit.py'

rule make_hub:
	'''
	creates a hub.txt file for viewing data on the UCSC genome browser
	'''
	input:
		mips_fasta=output+'/'+config['panel_name']+'_'+config['genome_name']+'_mips_filtered.psl'
	output:
		hub=output+'/'+config['panel_name']+'_'+config['genome_name']+'psl_hub.txt'
	script:
		'scripts/make_hub.py'
