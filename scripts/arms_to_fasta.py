arms_file=snakemake.input.arms_file
mips_fasta=open(snakemake.output.mips_fasta, 'w')
h_dict={}

def revcom(seq,nuc='DNA'):
	if nuc=='DNA':
		complement={'N':'N','n':'n','A':'T','C':'G','G':'C','T':'A','a':'t','t':'a','c':'g','g':'c', 'U':'A', 'u':'a', '-':'-'}
	else:
		complement={'N':'N','n':'n','A':'U','C':'G','G':'C','U':'A','a':'u','u':'a','c':'g','g':'c','-':'-'}
	return ''.join(reversed([complement[base] for base in seq]))

for line_number, line in enumerate(open(arms_file)):
	line=line.strip().split('\t')
	if line_number==0:
		for column_number, column_name in enumerate(line):
			h_dict[column_name]=column_number
		lig_c=h_dict['ligation_arm']
		ext_c=h_dict['extension_arm']
		gene_c=h_dict['gene_name']
		mip_family_c=h_dict['mip_family']
		mip_id_c=h_dict['mip_id']
	else:
		ext_seq=line[ext_c]
		lig_seq=revcom(line[lig_c])
		gene=line[gene_c]
		mip_family=line[mip_family_c]
		mip_id=line[mip_id_c]
		seq_name=f'>{gene}_{mip_family}_{mip_id}'
		catted_seq=ext_seq+lig_seq
		mips_fasta.write(seq_name+'\n')
		mips_fasta.write(catted_seq+'\n')
		
		
