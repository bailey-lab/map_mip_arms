'''
converts a bigbed file (in bigpsl format) into a track for viewing on the UCSC
genome browser.
'''
hub_name=snakemake.params.hub_name
email=snakemake.params.email
genome=snakemake.params.genome
big_bed=snakemake.params.big_bed
hub_file=open(snakemake.output.hub_file, 'w')

hub_file.write(f'hub {hub_name}\n')
hub_file.write(f'shortLabel {hub_name}\n')
hub_file.write(f'longLabel {hub_name}\n')
hub_file.write(f'useOneFile on\n')
hub_file.write(f'email {email}\n\n')
hub_file.write(f'genome {genome}\n\n')
hub_file.write(f'track {hub_name}\n')
hub_file.write(f'shortLabel {hub_name}\n')
hub_file.write(f'longLabel {hub_name}\n')
hub_file.write(f'visibility full\n')
hub_file.write(f'type bigPsl\n')
hub_file.write(f'bigDataUrl {big_bed}\n')
hub_file.write(f'color 40,40,40\n')