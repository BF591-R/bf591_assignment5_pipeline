trim_jar = '/projectnb/bubhub/conda_root/user_conda/jorofino/envs/r_for_bio/share/trimmomatic-0.39-2/trimmomatic.jar'
truseq_fa = '/projectnb/bubhub/conda_root/user_conda/jorofino/envs/r_for_bio/share/trimmomatic-0.39-2/adapters/TruSeq2-PE.fa'

import pandas
srrs = pandas.read_csv('sample_metadata.csv', usecols=['samplename'])['samplename'].tolist()

rule all:
	input:
		'samples_dl.done',
		'GRCm39.primary_assembly.genome.fa.gz',
		'gencode.vM27.primary_assembly.annotation.gtf.gz',
		'verse_counts.tsv'

rule wget_genome:
	output:
		'GRCm39.primary_assembly.genome.fa.gz'
	shell:
		'wget ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M27/GRCm39.primary_assembly.genome.fa.gz'

rule wget_gtf:
	output:
		'gencode.vM27.primary_assembly.annotation.gtf.gz'
	shell:
		'wget ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M27/gencode.vM27.primary_assembly.annotation.gtf.gz'

rule get_fastq_ebi_ena:
	params:
		url_list = 'fastq_ftp.txt'
	output:
		touch('samples_dl.done')
	shell:
		'wget -i {params.url_list}'

rule fastqc:
	input:
		fastq = '{srr}_{rp}.fastq.gz'
	output:
		fastqc = '{srr}_{rp, 1|2}_fastqc.zip'
	threads: 4
	shell:
		'fastqc -t {threads} {input.fastq}'

rule trimmomatic:
	input:
		r1 = '{srr}_1.fastq.gz',
		r2 = '{srr}_2.fastq.gz',
		trim_jar = trim_jar,
		truseq_fa = truseq_fa
	output:
		r1_p = '{srr}_1_P.fastq.gz',
		r1_u = '{srr}_1_U.fastq.gz',
		r2_p = '{srr}_2_P.fastq.gz',
		r2_u = '{srr}_2_U.fastq.gz',
	log: 
		'{srr}_trimlog.txt'
	threads: 8
	shell:
		'java -jar {input.trim_jar} PE '
		' -threads {threads} '
		' {input.r1} {input.r2} '
		' {output.r1_p} {output.r1_u} {output.r2_p} {output.r2_u} '
		' ILLUMINACLIP:{input.truseq_fa}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 2> {log} '

rule star:
	input:
		r1_p = '{srr}_1_P.fastq.gz',
		r2_p = '{srr}_2_P.fastq.gz',
		star_dir = 'GRCm39_primary_star_index'
	output:
		'{srr}.Aligned.out.bam'
	params:
		out_format = 'BAM Unsorted',
		prefix = '{srr}.',
	threads: 8
	shell:
		'STAR --genomeDir {input.star_dir} '
		' --runThreadN {threads} '
		' --readFilesIn {input.r1_p} {input.r2_p} '
		' --readFilesCommand zcat '
		' --outSAMtype {params.out_format} '
		' --outFileNamePrefix {params.prefix} '

rule verse:
	input:
		bams = '{srr}.Aligned.out.bam',
		gtf = 'gencode.vM27.primary_assembly.annotation.gtf'
	output:
		verse_out = '{srr}.exon.txt'
	params:
		prefix = '{srr}'
	threads: 4
	shell:
		'verse -T {threads} '
		' -a {input.gtf} '
		' -o {params.prefix} '
		' {input.bams} '

rule multiqc:
	input:
		expand('{srr}.exon.txt', srr=srrs)
	output:
		'multiqc_report.html'
	shell:
		'multiqc . -f'

rule concat_verse:
	input:
		meta = 'sample_metadata.csv'
	output:
		mat = 'verse_counts.tsv'
	run:
		import pandas
		meta = pandas.read_csv(input.meta)

		condition_d = dict(zip(meta['samplename'], meta['condition']))
		fns = meta['samplename'] + '.exon.txt' 


		dfs = []
		for fn in fns:
			condition = condition_d[fn.split('.')[0]]
			df = pandas.read_csv(fn, sep='\t', index_col='gene')
			df = df.rename(columns={'count':condition})
			dfs.append(df)
		
		pandas.concat(dfs, axis=1).to_csv(output.mat, sep='\t')
