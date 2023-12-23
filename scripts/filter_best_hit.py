mips_psl=snakemake.input.mips_psl
filtered_psl=open(snakemake.output.filtered_psl, 'w')
coords_file=open(snakemake.output.coords_file, 'w')

for line in open(mips_psl):
	split_line=line.strip().split('\t')
	if len(split_line)>4 and split_line[0].isdigit():
		diff=int(split_line[10])-int(split_line[0])
		if diff<1:
			filtered_psl.write(line)
	else:
		filtered_psl.write(line)
