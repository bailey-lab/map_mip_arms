configfile: 'unify_hub.yaml'

rule all:
	input:
		hub_file=config['output_folder']+'/mip_panel_hub.txt',
		snakemake=config['output_folder']+'/snakemake_params/unify_hub.smk'

rule copy_settings:
	input:
		snakemake='unify_hub.smk',
		yaml='unify_hub.yaml'
	output:
		snakemake=config['output_folder']+'/snakemake_params/unify_hub.smk',
		yaml=config['output_folder']+'/snakemake_params/unify_hub.yaml'
	shell:
		'''
		cp {input.snakemake} {output.snakemake}
		cp {input.yaml} {output.yaml}
		'''

rule create_hub:
	params:
		panels=config['panels'],
		hub_name=config['hub_name'],
		email=config['email'],
		genome=config['genome'],
		output_folder=config['output_folder']
	output:
		hub_file=config['output_folder']+'/mip_panel_hub.txt'
	script:
		'scripts/unify_hub.py'