import helper
import types

def sges():
	# assess the gene list
	for gene in args_.gene_list:
		if gene.lower() not in args_.target_genes:
			raise ValueError(f'Unrecognized target gene detected: {gene}')

	# update the arguments
	args_.gene_list = args_.target_genes if not args_.gene_list else [x.lower() for x in args_.gene_list]
	args_.target_regions = ' '.join([args_.gene_dict[x]['region'] for x in args_.gene_list])
	helper.create_dir(f'{args_.extended_prefix}.project')

	# write commands for each gene
	for target_gene in args_.gene_list:
		with open(f'{args_.extended_prefix}.project/{target_gene}.sh', 'w') as f:
			f.write(f'java -jar {args_.gatk} -T HaplotypeCaller -R {args_.assembly} -D {args_.dbsnp} --emitRefConfidence GVCF -U ALLOW_SEQ_DICT_INCOMPATIBILITY -I {args_.bam} -o {target_gene}.g.vcf -L {args_.gene_dict[target_gene]["region"]}\n')
			f.write(f'java -jar {args_.gatk} -T GenotypeGVCFs -R {args_.assembly} -D {args_.dbsnp} --variant {target_gene}.g.vcf -o {target_gene}.joint.vcf -L {args_.gene_dict[target_gene]["region"]}\n')
			f.write(f'java -jar {args_.gatk} -T VariantFiltration -R {args_.assembly} --filterExpression "QUAL <= 50.0" --filterName QUALFilter --variant {target_gene}.joint.vcf -o {target_gene}.joint.filtered.vcf -L {args_.gene_dict[target_gene]["region"]}\n')
			doc_region_list = [args_.gene_dict[target_gene]["region"], args_.control_region] if not args_.sampling else [args_.gene_dict[target_gene]["region"]] + helper.sample_regions(args_.control_region)
			doc_region_str = ' '.join([f'-L {x}' for x in helper.sort_regions(doc_region_list)])			
			f.write(f'java -jar {args_.gatk} -T DepthOfCoverage -R {args_.assembly} -I {args_.bam} --minMappingQuality {args_.mapping} --omitIntervalStatistics --omitPerSampleStats --omitLocusTable -U ALLOW_SEQ_DICT_INCOMPATIBILITY -o {target_gene}.gdf {doc_region_str}\n')
			f.write(f'python3 {args_.program_dir}/stargazer.py genotype -t {target_gene} -c {args_.control_gene} --vcf {target_gene}.joint.filtered.vcf --gdf {target_gene}.gdf -d wgs -o {target_gene}\n')

	# write report shell
	with open(f'{args_.extended_prefix}.project/report.sh', 'w') as f:
		shell = f'''head -n1 {args_.gene_list[0]}.stargazer-genotype.txt | awk '{{print "gene",$0}}' OFS="\\t" > merged.stargazer-genotype.txt

for gene in {' '.join(args_.gene_list)}
do
  tail -n+2 $gene.stargazer-genotype.txt | awk -v var="$gene" '{{print var,$0}}' OFS="\\t" >> merged.stargazer-genotype.txt
done

python3 {args_.program_dir}/stargazer.py report --genotype merged.stargazer-genotype.txt -o {args_.output_prefix}
'''
		f.write(shell)

	# write qsub shell
	with open(f'{args_.extended_prefix}.project/example-qsub.sh', 'w') as f:
		for target_gene in args_.gene_list:
			f.write(f'qsub -cwd -V -l mem_requested=30G -N {args_.output_prefix}-{target_gene} {target_gene}.sh\n')
		jid_list = [f'{args_.output_prefix}-{x}' for x in args_.gene_list]
		f.write(f'qsub -cwd -V -l mem_requested=30G -hold_jid {",".join(jid_list)} report.sh')

def read_manifest_file():
	if not any([args_.manifest, args_.sample_list, args_.bam_list]):
		raise TypeError(f"You must provide either argument 'minifest' or arguments 'sample_list' and 'bam_list'")
	if args_.manifest:
		if any([args_.sample_list, args_.bam_list]):
			raise TypeError(f"Argument 'manifest' cannot be used together with arguments 'sample_list' and 'bam_list'")
	else:
		if not all([args_.sample_list, args_.bam_list]):
			raise TypeError(f"Arguments 'sample_list' and 'bam_list' require each other")
		if len(args_.sample_list) != len(args_.bam_list):
			raise TypeError(f"Arguments 'sample_list' and 'bam_list' must be of same length")

	args_.batch_list = []
	args_.sample_lists = []

	for i in range(args_.batch):
		args_.batch_list.append('b{}'.format(str(i + 1).zfill(len(str(args_.batch)))))
	sample_list = []
	
	if args_.manifest:
		with open(args_.manifest) as f:
			header = next(f).strip().split('\t')
			s = header.index('sample_id')
			b = header.index('bam')
			for line in f:
				fields = line.strip().split('\t')
				sample = types.SimpleNamespace()
				sample.sample_id = fields[s]
				sample.bam = fields[b]
				sample_list.append(sample)
	else:
		for i in range(len(args_.sample_list)):
			sample = types.SimpleNamespace()
			sample.sample_id = args_.sample_list[i]
			sample.bam = args_.bam_list[i]
			sample_list.append(sample)
			
	s = len(sample_list)
	b = args_.batch
	args_.sample_lists = [sample_list[(x * s // b):((x + 1) * s // b)] for x in range(b)]
	with open(args_.log, 'a') as f:
		f.write(f'\n{args_.line_break}\nGathering project summary...\n\nStatus: Completed\n\nSamples total: {s}\nBatches total: {b}\n\n{args_.line_break}')

	helper.create_dir(f'{args_.extended_prefix}.project')
	for b in args_.batch_list:
		helper.create_dir(f'{args_.extended_prefix}.project/{b}')

def write_rbind_shell():
	with open(f'{args_.extended_prefix}.project/rbind.sh', 'w') as f:
		shell = f'''head -n1 {args_.batch_list[0]}.stargazer-genotype.txt | awk '{{print "batch",$0}}' OFS="\\t" > merged.stargazer-genotype.txt

for b in {' '.join(args_.batch_list)}
do
  tail -n+2 $b.stargazer-genotype.txt | awk -v var="$b" '{{print var,$0}}' OFS="\\t" >> merged.stargazer-genotype.txt
done

python3 {args_.program_dir}/stargazer.py view summary --genotype merged.stargazer-genotype.txt -o merged
'''
		f.write(shell)

def write_qsub_shell(pars):
	with open(f'{args_.extended_prefix}.project/example-qsub.sh', 'w') as f:
		for i, b in enumerate(args_.batch_list):
			f.write(f'qsub -cwd -V {pars} -e {b} -o {b} -N {b}-doc.{args_.output_prefix} {b}/doc.sh\n')
			for sample in args_.sample_lists[i]:
				f.write(f'qsub -cwd -V {pars} -e {b} -o {b} -N {b}-hc.{args_.output_prefix} {b}/{sample.sample_id}.hc.sh\n')
			f.write(f'qsub -cwd -V {pars} -e {b} -o {b} -hold_jid {b}-hc.{args_.output_prefix} -N {b}-posthc.{args_.output_prefix} {b}/posthc.sh\n')
			f.write(f'qsub -cwd -V {pars} -e {b} -o {b} -hold_jid {b}-posthc.{args_.output_prefix},{b}-doc.{args_.output_prefix} -N {b}-rs.{args_.output_prefix} {b}/rs.sh\n')
		_ = ','.join([f'{x}-rs.{args_.output_prefix}' for x in args_.batch_list])
		f.write(f'qsub -cwd -V {pars} -hold_jid {_} -N rbind.{args_.output_prefix} rbind.sh\n')

def sgep():
	read_manifest_file()

	doc_region_list = [args_.target_region, args_.control_region] if not args_.sampling else [args_.target_region] + helper.sample_regions(args_.control_region)
	doc_region_str = ' '.join([f'-L {x}' for x in helper.sort_regions(doc_region_list)])

	for i, b in enumerate(args_.batch_list):
		with open(f'{args_.extended_prefix}.project/{b}/bam.list', 'w') as f:
			for sample in args_.sample_lists[i]:
				f.write(sample.bam + '\n')

		with open(f'{args_.extended_prefix}.project/{b}/doc.sh', 'w') as f:
			f.write(f'java -jar {args_.gatk} -T DepthOfCoverage -R {args_.assembly} -I {b}/bam.list --minMappingQuality {args_.mapping} --omitIntervalStatistics --omitPerSampleStats --omitLocusTable -U ALLOW_SEQ_DICT_INCOMPATIBILITY -o {b}.gdf {doc_region_str}\n')

		for sample in args_.sample_lists[i]:
			with open(f'{args_.extended_prefix}.project/{b}/{sample.sample_id}.hc.sh', 'w') as f:
				f.write(f'java -jar {args_.gatk} -T HaplotypeCaller -R {args_.assembly} -D {args_.dbsnp} --emitRefConfidence GVCF -U ALLOW_SEQ_DICT_INCOMPATIBILITY -I {sample.bam} -o {b}/{sample.sample_id}.g.vcf -L {args_.target_region}\n')
		sample_id_str = ' '.join([f'--variant {b}/{x.sample_id}.g.vcf' for x in args_.sample_lists[i]])

		with open(f'{args_.extended_prefix}.project/{b}/posthc.sh', 'w') as f:
			f.write(f'java -jar {args_.gatk} -T GenotypeGVCFs -R {args_.assembly} -D {args_.dbsnp} -o {b}/{args_.output_prefix}.joint.vcf -L {args_.target_region} {sample_id_str}\n')
			f.write(f'java -jar {args_.gatk} -T VariantFiltration -R {args_.assembly} --filterExpression "QUAL <= 50.0" --filterName QUALFilter --variant {b}/{args_.output_prefix}.joint.vcf -o {b}.joint.filtered.vcf -L {args_.target_region}\n')

		with open(f'{args_.extended_prefix}.project/{b}/rs.sh', 'w') as f:
			f.write(f'python3 {args_.program_dir}/stargazer.py genotype -t {args_.target_gene} -c {args_.control_gene} --vcf {b}.joint.filtered.vcf --gdf {b}.gdf -d {args_.data_type} -o {b}\n')

	write_rbind_shell()
	write_qsub_shell('-l mem_requested=30G')

def sgea():
	read_manifest_file()

	for i, b in enumerate(args_.batch_list):
		with open(f'{args_.extended_prefix}.project/{b}/bam.list', 'w') as f:
			for sample in args_.sample_lists[i]:
				f.write(sample.bam + '\n')
				
		doc_region_list = [args_.target_region, args_.control_region] if not args_.sampling else [args_.target_region] + helper.sample_regions(args_.control_region)
		
		with open(f'{args_.extended_prefix}.project/{b}/doc.sh', 'w') as f:
			for region in helper.sort_regions(doc_region_list):
				f.write(f'samtools depth -a -Q {args_.mapping} -r chr{region} -f {b}/bam.list >> {b}/{b}.sdf\n')
			f.write(f'python3 {args_.program_dir}/stargazer.py setup sdf2gdf --sdf {b}/{b}.sdf -o {b} --sample_list {" ".join([x.sample_id for x in args_.sample_lists[i]])}\n')

		for sample in args_.sample_lists[i]:
			with open(f'{args_.extended_prefix}.project/{b}/{sample.sample_id}.hc.sh', 'w') as f:
				f.write(f'gatk HaplotypeCaller -R {args_.assembly} -D {args_.dbsnp} --emit-ref-confidence GVCF -I {sample.bam} -O {b}/{sample.sample_id}.g.vcf -L chr{args_.target_region}\n')		

		with open(f'{args_.extended_prefix}.project/{b}/posthc.sh', 'w') as f:
			variant_list_str = ' '.join([f'--variant {b}/{x.sample_id}.g.vcf' for x in args_.sample_lists[i]])
			f.write(f'gatk CombineGVCFs -R {args_.assembly} -D {args_.dbsnp} -O {b}/{b}.g.vcf -L chr{args_.target_region} {variant_list_str}\n')
			f.write(f'gatk GenotypeGVCFs -R {args_.assembly} -D {args_.dbsnp} -O {b}/{b}.joint.vcf -L chr{args_.target_region} --variant {b}/{b}.g.vcf\n')
			f.write(f'gatk VariantFiltration -R {args_.assembly} --filter-expression "QUAL <= 50.0" --filter-name QUALFilter --variant {b}/{b}.joint.vcf -O {b}.joint.filtered.vcf -L chr{args_.target_region}\n')

		with open(f'{args_.extended_prefix}.project/{b}/rs.sh', 'w') as f:
			f.write(f'python3 {args_.program_dir}/stargazer.py genotype -t {args_.target_gene} -c {args_.control_gene} --vcf {b}.joint.filtered.vcf --gdf {b}.gdf -d {args_.data_type} -o {b}\n')

	write_rbind_shell()
	write_qsub_shell('-q biall.q -S /bin/bash')

def fastq2bam():
	if not args_.manifest:
		raise TypeError('Argument "--manifest" not given')
	if not args_.assembly:
		raise TypeError('Argument "--assembly" not given')
	if not args_.variant_list:
		raise TypeError('Argument "--variant_list" not given')
	if not args_.bed:
		raise TypeError('Argument "--bed" not given')
	
	sample_dict = {}

	with open(args_.manifest) as f:
		header = next(f).strip().split('\t')
		for line in f:
			fields = line.strip().split('\t')
			sample_id = fields[header.index('sample_id')]
			fastq = fields[header.index('fastq')]
			if sample_id not in sample_dict:
				sample_dict[sample_id] = []
			sample_dict[sample_id].append(fastq)
			if len(sample_dict[sample_id]) > 2:
				raise ValueError('Check your manifest file')
	
	shell_dir = f'{args_.extended_prefix}.project/shell'
	bam_dir = f'{args_.extended_prefix}.project/bam'
	eno_dir = f'{args_.extended_prefix}.project/eno'
	temp_dir = f'{args_.extended_prefix}.project/temp'
	stat_dir = f'{args_.extended_prefix}.project/stat'
	
	helper.create_dir(f'{args_.extended_prefix}.project')
	helper.create_dir(shell_dir)
	helper.create_dir(bam_dir)
	helper.create_dir(eno_dir)
	helper.create_dir(temp_dir)
	helper.create_dir(stat_dir)
	
	for sample_id in sample_dict:
		with open(f'{shell_dir}/{sample_id}.s1.sh', 'w') as f:
			f.write('#!/bin/bash\n')
			f.write('\n')
			f.write(f'name={sample_id}\n')
			f.write(f'read1={sample_dict[sample_id][0]}\n')
			f.write(f'read2={sample_dict[sample_id][1]}\n')
			f.write(f'thread={args_.thread}\n')
			f.write(f'assembly={args_.assembly}\n')
			f.write(f'platform={args_.platform}\n')
			f.write(f'library={args_.output_prefix}\n')
			f.write(f'bam1={args_.output_dir}/{temp_dir}/$name.sorted.bam\n')
			f.write('\n')
			f.write('# Get read group information\n')
			f.write('first=`zcat $read1 | head -1`\n')
			f.write('''flowcell=`echo "$first" | awk -F " " '{print $1}' | awk -F ":" '{print $3}'`\n''')
			f.write('''barcode=`echo "$first" | awk -F " " '{print $2}' | awk -F ":" '{print $4}'`\n''')
			f.write(f'group="@RG\\tID:$flowcell\\tPU:$flowcell.$barcode\\tSM:$name\\tPL:$platform\\tLB:$library"\n')
			f.write('\n')
			f.write('# Align and sort seuqnece reads\n')
			f.write('# Note that this step will also assign read group to mapped reads\n')
			f.write('bwa mem -M -R $group -t $thread $assembly $read1 $read2 | samtools sort -@ $thread -o $bam1 -\n')

		with open(f'{shell_dir}/{sample_id}.s2.sh', 'w') as f:
			f.write('#!/bin/bash\n')
			f.write('\n')
			f.write(f'name={sample_id}\n')
			f.write(f'length={args_.length}\n')
			f.write(f'assembly={args_.assembly}\n')
			if args_.bed:
				f.write(f'bed={args_.bed}\n')
			f.write(f'bam1={args_.output_dir}/{temp_dir}/$name.sorted.bam\n')
			f.write(f'bam2={args_.output_dir}/{temp_dir}/$name.sorted.markeddups.bam\n')
			f.write(f'bam3={args_.output_dir}/{bam_dir}/$name.sorted.markeddups.recal.bam\n')
			f.write(f'table={args_.output_dir}/{temp_dir}/$name.table\n')
			f.write(f'metrics={args_.output_dir}/{temp_dir}/$name.metrics\n')
			f.write(f'stat={args_.output_dir}/{stat_dir}/$name.txt\n')
			
			n = 0
			site_list = []
			for variant in args_.variant_list:
				n += 1
				f.write(f'site{n}={variant}\n')
				site_list.append(f'$site{n}')

			f.write('\n')
			f.write('# Mark duplicate reads\n')
			f.write(f'gatk MarkDuplicates -I $bam1 -M $metrics -O $bam2 --ASSUME_SORT_ORDER coordinate\n')
			f.write('samtools index $bam2\n')
			f.write('\n')
			f.write('# Build BQSR model\n')
			f.write(f'gatk BaseRecalibrator -I $bam2 -O $table -R $assembly --intervals $bed {" ".join(["--known-sites "+ x for x in site_list])}\n')
			f.write('\n')
			f.write('# Apply BQSR model\n')
			f.write('# Note that this step will remove mapped reads outside of targeted regions\n')
			f.write(f'gatk ApplyBQSR -I $bam2 -O $bam3 --intervals $bed -bqsr $table\n')
			f.write('\n')
			f.write('# Get various statistics\n')
			f.write('flagstat1=`samtools flagstat $bam2`\n')
			f.write('flagstat2=`samtools flagstat $bam3`\n')
			f.write('''total_n=`echo "$flagstat1" | awk 'NR == 1 {print $1}'`\n''')
			f.write('''mapped_n=`echo "$flagstat1" | awk 'NR == 5 {print $1}'`\n''')
			f.write('''mapped_p=`echo "$flagstat1" | awk 'NR == 5 {print $5}' | sed "s/(//g" | sed "s/%//g"`\n''')
			f.write('''unique_n=`samtools view -q 1 -c $bam2`\n''')
			f.write('''unique_p=`bc <<< "scale = 2; $unique_n * 100 / $total_n"`\n''')
			f.write('''duplicate_n=`echo "$flagstat1" | awk 'NR == 4 {print $1}'`\n''')
			f.write('''duplicate_p=`bc <<< "scale = 2; $duplicate_n * 100 / $total_n"`\n''')
			f.write('''total_b=`echo "$(($total_n * $length))"`\n''')
			f.write('''panel_b=`awk -F '\\t' 'BEGIN {sum = 0} {sum += $3 - $2} END {print sum}' $bed`\n''')
			f.write('''total_x=$(($total_b / $panel_b))\n''')
			f.write('''target_n=`echo "$flagstat2" | awk 'NR == 1 {print $1}'`\n''')
			f.write('''target_p=`bc <<< "scale = 2; $target_n * 100 / $total_n"`\n''')
			f.write('''target_b=`echo "$(($target_n * $length))"`\n''')
			f.write('''target_x=`echo "$(($target_b / $panel_b))"`\n''')
			f.write('''mean_x=`samtools depth $bam3 -b $bed | awk '{sum += $3} END {print sum / NR}' | sed 's/\.[^\.]*$//'`\n''')
			f.write('''x0001_p=`printf "%.2f\\n" <<< echo "$(samtools depth -b $bed $bam3 | awk '{if ($3 >= 1) c++} END {print c / NR * 100}')"`\n''')
			f.write('''x0010_p=`printf "%.2f\\n" <<< echo "$(samtools depth -b $bed $bam3 | awk '{if ($3 >= 10) c++} END {print c / NR * 100}')"`\n''')
			f.write('''x0020_p=`printf "%.2f\\n" <<< echo "$(samtools depth -b $bed $bam3 | awk '{if ($3 >= 20) c++} END {print c / NR * 100}')"`\n''')
			f.write('''x0050_p=`printf "%.2f\\n" <<< echo "$(samtools depth -b $bed $bam3 | awk '{if ($3 >= 50) c++} END {print c / NR * 100}')"`\n''')
			f.write('''x0100_p=`printf "%.2f\\n" <<< echo "$(samtools depth -b $bed $bam3 | awk '{if ($3 >= 100) c++} END {print c / NR * 100}')"`\n''')
			f.write('''x0200_p=`printf "%.2f\\n" <<< echo "$(samtools depth -b $bed $bam3 | awk '{if ($3 >= 200) c++} END {print c / NR * 100}')"`\n''')
			f.write('''x0300_p=`printf "%.2f\\n" <<< echo "$(samtools depth -b $bed $bam3 | awk '{if ($3 >= 300) c++} END {print c / NR * 100}')"`\n''')
			f.write('''x0400_p=`printf "%.2f\\n" <<< echo "$(samtools depth -b $bed $bam3 | awk '{if ($3 >= 400) c++} END {print c / NR * 100}')"`\n''')
			f.write('''x0500_p=`printf "%.2f\\n" <<< echo "$(samtools depth -b $bed $bam3 | awk '{if ($3 >= 500) c++} END {print c / NR * 100}')"`\n''')
			f.write('''x1000_p=`printf "%.2f\\n" <<< echo "$(samtools depth -b $bed $bam3 | awk '{if ($3 >= 1000) c++} END {print c / NR * 100}')"`\n''')
			f.write('\n')
			f.write('# Report various statistics\n')
			f.write('''echo "total_n $total_n" > $stat\n''')
			f.write('''echo "mapped_n $mapped_n" >> $stat\n''')
			f.write('''echo "mapped_p $mapped_p" >> $stat\n''')
			f.write('''echo "unique_n $unique_n" >> $stat\n''')
			f.write('''echo "unique_p $unique_p" >> $stat\n''')
			f.write('''echo "duplicate_n $duplicate_n" >> $stat\n''')
			f.write('''echo "duplicate_p $duplicate_p" >> $stat\n''')
			f.write('''echo "total_b $total_b" >> $stat\n''')
			f.write('''echo "panel_b $panel_b" >> $stat\n''')
			f.write('''echo "total_x $total_x" >> $stat\n''')
			f.write('''echo "target_n $target_n" >> $stat\n''')
			f.write('''echo "target_p $target_p" >> $stat\n''')
			f.write('''echo "target_b $target_b" >> $stat\n''')
			f.write('''echo "target_x $target_x" >> $stat\n''')
			f.write('''echo "mean_x $mean_x" >> $stat\n''')
			f.write('''echo "x0001_p $x0001_p" >> $stat\n''')
			f.write('''echo "x0010_p $x0010_p" >> $stat\n''')
			f.write('''echo "x0020_p $x0020_p" >> $stat\n''')
			f.write('''echo "x0050_p $x0020_p" >> $stat\n''')
			f.write('''echo "x0100_p $x0100_p" >> $stat\n''')
			f.write('''echo "x0200_p $x0200_p" >> $stat\n''')
			f.write('''echo "x0300_p $x0300_p" >> $stat\n''')
			f.write('''echo "x0400_p $x0400_p" >> $stat\n''')
			f.write('''echo "x0500_p $x0500_p" >> $stat\n''')
			f.write('''echo "x1000_p $x1000_p" >> $stat\n''')

	with open(f'{args_.extended_prefix}.project/example-qsub.sh', 'w') as f:
		f.write(f'eno={args_.output_dir}/{eno_dir}\n')
		f.write(f'shell={args_.output_dir}/{shell_dir}\n')
		f.write('\n')
		for sample_id in sample_dict:
			f.write(f'qsub -V -q biall.q -S /bin/bash -e $eno -o $eno -N {sample_id}.s1 -pe pePAC {args_.thread}  $shell/{sample_id}.s1.sh\n')
			f.write(f'qsub -V -q biall.q -S /bin/bash -e $eno -o $eno -N {sample_id}.s2 -hold_jid {sample_id}.s1 $shell/{sample_id}.s2.sh\n')

SUBTOOLS = {'sges': sges, 'sgep': sgep, 'sgea': sgea, 'fastq2bam': fastq2bam}
DESCRIPTION = f'''tool description:
  provide end-to-end solutions for genotyping pipeline

getting help:
  stargazer.py pipeline -h

available subtools:
  sges
    run per-sample genotyping pipeline (i.e. single target sample, multiple target genes) with SGE and GATK3
    
  sgep
    run per-project genotyping pipeline (i.e. multiple target samples, single target gene) with SGE and GATK3
    
  sgea
    run per-project genotyping pipeline (i.e. multiple target samples, single target gene) with SGE and Anaconda (more specifically, with GATK4 and SAMtools)

  fastq2bam
    run fastq2bam pipeline with SGE and Anaconda (more specifically, with BWA-MEM, GATK4, and SAMtools)

main usages:
  run the sges subtool
    stargazer.py pipeline sges -o OUTPUT_PREFIX --assembly ASSEMBLY -c CONTROL_GENE --dbsnp DBSNP --gatk GATK --bam BAM [--gene_list [GENE [GENE ...]]]
  
  run the sgep subtool with manifest file
    stargazer.py pipeline sgep -o OUTPUT_PREFIX --assembly ASSEMBLY -c CONTROL_GENE --dbsnp DBSNP -d DATA_TYPE --gatk GATK --manifest MANIFEST -t TARGET_GENE [--batch BATCH]
  
  run the sgea subtool with manifest file
    stargazer.py pipeline sgea -o OUTPUT_PREFIX --assembly ASSEMBLY -c CONTROL_GENE --dbsnp DBSNP -d DATA_TYPE --manifest MANIFEST -t TARGET_GENE [--batch BATCH]

  run the fastq2bam subtool with manifest file
    stargazer.py pipeline fastq2bam -o OUTPUT_PREFIX --assembly ASSEMBLY --bed BED --manifest MANIFEST [--thread THREAD] --variant_list [VARIANT [VARIANT ...]]

other usages:
  run the sgep subtool without manifest file
    stargazer.py pipeline sgep -o OUTPUT_PREFIX --assembly ASSEMBLY -c CONTROL_GENE --dbsnp DBSNP -d DATA_TYPE --gatk GATK --sample_list [SAMPLE [SAMPLE ...]] --bam_list [BAM [BAM ...]] -t TARGET_GENE [--batch BATCH]
  
  run the sgea subtool without manifest file
    stargazer.py pipeline sgea -o OUTPUT_PREFIX --assembly ASSEMBLY -c CONTROL_GENE --dbsnp DBSNP -d DATA_TYPE --sample_list [SAMPLE [SAMPLE ...]] --bam_list [BAM [BAM ...]] -t TARGET_GENE [--batch BATCH]
'''

def run(args):
	global args_
	args_ = args

	SUBTOOLS[args_.subtool]()
