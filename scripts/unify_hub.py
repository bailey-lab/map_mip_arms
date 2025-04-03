'''
combines all MIP panels together into a single hub with multiple tracks.
'''
import os
import subprocess

panels=snakemake.params.panels
hub_name=snakemake.params.hub_name
email=snakemake.params.email
genome=snakemake.params.genome
output_folder=snakemake.params.output_folder
output_hub_file=open(snakemake.output.hub_file, 'w')

output_hub_file.write(f'hub {hub_name}\n')
output_hub_file.write(f'shortLabel {hub_name}\n')
output_hub_file.write(f'longLabel {hub_name}\n')
output_hub_file.write(f'useOneFile on\n')
output_hub_file.write(f'email {email}\n\n')
output_hub_file.write(f'genome {genome}\n\n')

for panel in panels:
	panel_files=os.listdir(panel)
	for panel_file in panel_files:
		if panel_file.endswith('hub.txt'):
			input_hub_file=panel_file
		elif panel_file.endswith('bigPsl.bb'):
			psl_file=panel_file
	subprocess.call(['cp', panel+'/'+psl_file, output_folder])
	for line_number, line in enumerate(open(panel+'/'+input_hub_file)):
		if line_number>7:
			output_hub_file.write(line)
	output_hub_file.write('\n')