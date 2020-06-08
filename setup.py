import helper
import statistics
import os

def sdf2gdf():
	header = ['Locus', 'Total_Depth', 'Average_Depth_sample'] + ['Depth_for_' + x for x in args_.sample_list]

	# peek into the first line
	with open(args_.sdf) as f:
		fields = next(f).strip().split('\t')
		if len(fields) != len(header) - 1:
			raise ValueError('Input SDF and sample list different in length')

	with open(args_.output_prefix + '.gdf', 'w') as f1:
		f1.write('\t'.join(header) + '\n')
		with open(args_.sdf) as f2:
			for line in f2:
				fields = line.strip().split('\t')
				chr = fields[0].replace('chr', '')
				pos = fields[1]
				locus = chr + ':' + pos
				data = [int(x) for x in fields[2:]]
				total = sum(data)
				avg = round(statistics.mean(data), 2)
				new_fields = [locus, total, avg] + data
				f1.write('\t'.join([str(x) for x in new_fields]) + '\n')

def define():
	with open(args_.table) as f:
		s = 14 # first star allele index
		header = next(f).strip().split('\t')
		dat = [[] for x in header[s:]]
		for line in f:
			fields = line.strip().split('\t')
			pos = fields[5]
			var = fields[8]
			wt = fields[9]
			for i in range(s, len(fields)):
				if fields[i] == '1':
					dat[i - s].append('{}:{}>{}'.format(pos, wt, var))

	with open(args_.output_prefix + '.stargazer-view.txt', 'w') as f:
		for i in range(len(dat)):
			f.write(header[s + i] + '\t' + ','.join(dat[i]) + '\n')

def slice():
	sliced_vcf = helper.read_vcf_region(args_.vcf, args_.region)
	helper.write_vcf(sliced_vcf, args_.output_prefix + '.vcf')

def check():
	with open(args_.log, 'a') as f:
		f.write(f'\n{args_.line_break}\n')
		f.write('Step 1/2: Checking star_table.txt...\n\n')

	star_data = []
	snp_data = []
	gene_list = []
	
	with open(f'{args_.program_dir}/star_table.txt') as f:
		star_header = next(f).strip().split('\t')
		for line in f:
			fields = line.strip().split('\t')
			gene = fields[star_header.index('gene')]
			if gene not in gene_list:
				gene_list.append(gene)			
			star_data.append(fields)
			
	with open(f'{args_.program_dir}/snp_table.txt') as f:
		snp_header = next(f).strip().split('\t')
		for line in f:
			fields = line.strip().split('\t')
			snp_data.append(fields)
	
	snp_dict = {}
	star_dict = {}
	
	for gene in gene_list:
		snp_dict[gene] = [x for x in snp_data if x[snp_header.index('gene')] == gene]
		star_dict[gene] = [x for x in star_data if x[star_header.index('gene')] == gene]
		
	def check_allele(x, gene, has_cadd30, type, has_lof):
		snp_list = [x[snp_header.index('pos')] + ':' + x[snp_header.index('wt')] + '>' + x[snp_header.index('var')] for x in snp_dict[gene]]
		cadd_dict = dict(zip(snp_list, [float(x[snp_header.index('cadd')]) for x in snp_dict[gene]]))
		so_dict = dict(zip(snp_list, [x[snp_header.index('so')] for x in snp_dict[gene]]))
		effect_dict = dict(zip(snp_list, [x[snp_header.index('effect')] for x in snp_dict[gene]]))
		lof_dict = dict(zip(snp_list, [x[snp_header.index('lof')] for x in snp_dict[gene]]))
		seen_list = []	
		if x == '.' or x == 'ref':
			return
			
		cadd30_seen = False
			
		for snp in x.split(','):
			pos = int(snp.split(':')[0])
			wt = snp.split(':')[1].split('>')[0]
			var = snp.split(':')[1].split('>')[1]
			
			if snp not in snp_list:
				raise ValueError(f'Unrecognized allele definition: {gene.upper()}{name} (#{number}) has {snp}')
			if wt == var:
				raise ValueError(f'Incorrect allele definition: {gene.upper()}{name} (#{number}) has {pos}-{wt}-{var}')
			if seen_list and seen_list[-1] > pos:
				raise ValueError(f'Allele definition is not coordinate sorted: {gene.upper()}{name} (#{number}) has {x}')
			seen_list.append(pos)
			
			if type == 'core':
				if cadd_dict[snp] >= 30:
					cadd30_seen = True
					if has_cadd30 == 'no':
						raise ValueError(f'Incorrect CADD information: {gene.upper()}{name} (#{number}) is marked as not having CADD30, but it has {snp} with CADD >= 30 ({cadd_dict[snp]})')
				
				if has_lof == 'no' and lof_dict[snp] == 'yes':
					raise ValueError(f'Incorrect LoF information: {gene.upper()}{name} (#{number}) is marked as not having LoF, but it has {snp} with LoF (Effect={effect_dict[snp]}; SO={so_dict[snp]})')
			
			
		if type == 'core' and not cadd30_seen and has_cadd30 == 'yes':
			raise ValueError(f'Incorrect CADD information: {gene.upper()}{name} (#{number}) is marked as having CADD30, but no such variant was found ({x}; {",".join({str(cadd_dict[y]) for y in x.split(",")})})')

	for fields in star_data:
		gene = fields[star_header.index('gene')]
		name = fields[star_header.index('name')]
		core = fields[star_header.index('core')]
		tag = fields[star_header.index('tag')]
		number = int(fields[star_header.index('number')])
		has_cadd30 = fields[star_header.index('has_cadd30')]
		has_lof = fields[star_header.index('has_lof')]
		score = fields[star_header.index('score')]
		if number == 1:
			count = 0
		count += 1
		if number != count:
			raise ValueError(f'Unmatched allele count: {gene.upper()}{name} (#{number}) should count {count}')
		check_allele(core, gene, has_cadd30, 'core', has_lof)
		check_allele(tag, gene, has_cadd30, 'tag', has_lof)
		if has_cadd30 == 'yes' and score == 'unknown':
			raise ValueError(f'Conflicting data: {gene.upper()}{name} (#{number}) is marked as having CADD30, but it has unknown score')
		if has_lof == 'yes' and score == 'unknown':
			raise ValueError(f'Conflicting data: {gene.upper()}{name} (#{number}) is marked as having LoF, but it has unknown score')

	with open(args_.log, 'a') as f:
		f.write(f'Status: Completed\n\n')
		f.write('Allele counts: Pass\n')
		f.write('Allele definition: Pass\n')
		f.write('CADD30 information: Pass\n')
		f.write('LoF information: Pass\n')
		f.write(f'\n{args_.line_break}\n')

	with open(args_.log, 'a') as f:
		f.write('Step 2/2: Checking snp_table.txt...\n\n')

	for fields in snp_data:			
		gene = fields[snp_header.index('gene')]
		pos = fields[snp_header.index('pos')]
		hg = fields[snp_header.index('hg')]
		var = fields[snp_header.index('var')]
		wt = fields[snp_header.index('wt')]
		rev = fields[snp_header.index('rev')] == 'yes'
		impact = fields[snp_header.index('impact')]
		cadd = float(fields[snp_header.index('cadd')])
		causal = fields[snp_header.index('causal')]
		
		name = f'{gene.upper()}-{pos}-{hg}-{var}-{wt}'
		number = int(fields[1])
		current_pos = int(pos)
		if number == 1:
			previous_pos = current_pos
			previous_var = var
			count = 0
		count += 1
		if current_pos < previous_pos:
			raise ValueError(f'Variants are not coordinate sorted: {name} (#{number}) is after {previous_pos}')
		if current_pos == previous_pos and var < previous_var:
			raise ValueError(f'Variants are not alphabetically sorted: {name} is after {previous_var}')
		previous_pos = current_pos
		previous_var = var
		if number != count:
			raise ValueError(f'Unmatched variant count: {name} (#{number}) should count {count}')
		if rev and (hg != var or wt == hg or wt == var):
			raise ValueError(f'Conflicting data: {name} is marked as revertant')
		if cadd >= 30 and impact != 'high_impact':
			raise ValueError(f'Conflicting data: {name} has CADD >= 30 and {impact}')
		
		# check associated phenotype
		phenotype = ''
		for star in star_dict[gene]:
			if ',' in star[star_header.index('core')]:
				continue
			if star[star_header.index('sv')] != '.':
				continue
			if f'{pos}:{wt}>{var}' in star[star_header.index('core')].split(','):
				phenotype = star[star_header.index('phenotype')]
		if not phenotype:
			phenotype = '.'		
		if causal != phenotype:
			raise ValueError(f'Incorrect phenotype information: {name} ({",".join(causal)}) should have {phenotype}')
		if causal not in ['.', 'unknown_function', 'normal_function', 'IV/Normal'] and impact != 'high_impact':
			raise ValueError(f'Incorrect impact information: {name} ({impact}; {causal}) should have high_impact')

	with open(args_.log, 'a') as f:
		f.write(f'Status: Completed\n\n')
		f.write('Variant counts: Pass\n')
		f.write('Coordinate sorted: Pass\n')
		f.write('Alphabetically sorted: Pass\n')
		f.write('Reverting variants: Pass\n')
		f.write('Variant impact: Pass\n')
		f.write('Variant phenotype: Pass\n')
		f.write(f'\n{args_.line_break}')

def merge():
	if not args_.vcf_dir:
		raise ValueError('Required argument "--vcf_dir" not found')

	vcfs = []
	region = None

	if args_.target_gene and args_.region:
		raise ValueError('Arguments "--target_gene" and "--region" cannot be used together')
	if args_.target_gene:
		region = helper.get_region(args_.target_gene)
	if args_.region:
		region = args_.region

	for r, d, f in os.walk(args_.vcf_dir):
		for file in f:
			if file.endswith('vcf'):
				vcf_path = os.path.join(r, file)
				if region:
					vcf = helper.read_vcf_region(vcf_path, region)
				else:
					vcf = helper.read_vcf_simple(vcf_path)
				vcfs.append(vcf)
		
	merged_vcf = vcfs[0]
	size = len(merged_vcf.data)

	for vcf in vcfs[1:]:
		sample_id = vcf.header[9]
		merged_vcf.header.append(sample_id)
		for i in range(size):
			data = vcf.data[i][9]
			merged_vcf.data[i].append(data)

	helper.write_vcf(merged_vcf, f'{args_.output_prefix}.vcf')

SUBTOOLS = {'sdf2gdf': sdf2gdf, 'define': define, 'slice': slice, 'check': check, 'merge': merge}

DESCRIPTION = f'''tool description:
  create various files necessary for running Stargazer

getting help:
  stargazer.py setup -h

available subtools:
  sdf2gdf
    create gdf file from sdf file

  define
    define star alleles using variants from the SNP table
    
  slice
    create sliced VCF file

  check
    check the SNP and star allele tables

  merge
    create new VCF file by merging multiple VCF files

main usages:
  run the sdf2gdf subtool
    stargazer.py setup sdf2gdf -o OUTPUT_PREFIX --sample_list [SAMPLE [SAMPLE ...]] --sdf SDF

  run the define subtool
    stargazer.py setup define -o OUTPUT_PREFIX --table TABLE

  run the slice subtool
    stargazer.py setup slice -o OUTPUT_PREFIX --region REGION --VCF VCF

  run the check subtool
    stargazer.py setup check -o OUTPUT_PREFIX

  run the merge subtool
    stargazer.py setup merge -o OUTPUT_PREFIX --vcf_dir VCF_DIR [--region REGION] [-t TARGET_GENE]
'''

def run(args):
	global args_
	args_ = args

	SUBTOOLS[args_.subtool]()