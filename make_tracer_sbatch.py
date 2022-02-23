#! /usr/bin/python3

import argparse, os

def find_fastq(path):
	for name in os.listdir(path):
		path2 = os.path.join(path, name)
		if os.path.isdir(path2):
			for fastq in find_fastq(path2): yield fastq
		elif any(name.endswith(ending) for ending in ('fastq.gz','fastq','fq.gz','fq')):
			yield path2


# command: python3 make_tracer_sbatch.py trimmed_fastq/ tracer-sbatch/

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('folder_in')
	parser.add_argument('folder_out')
	o = parser.parse_args()
	
	o.folder_in = os.path.realpath(o.folder_in)
	o.folder_out = os.path.realpath(o.folder_out)
	
	for sample in os.listdir(o.folder_in):
		file_out = os.path.join(o.folder_out, sample+'_tracer.sh')
		fastqs = list(find_fastq(os.path.join(o.folder_in, sample)))
		if not fastqs:
			print(os.path.join(o.folder_in, sample))
			raise Exception
		with open(file_out, 'w') as outfh:
			print('''#!/bin/bash -l
#SBATCH -A sens2018579
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1:00:00''', file=outfh)
			print('#SBATCH -J tracer_'+sample, file=outfh)
			print('''
module load bioinfo-tools
module load bowtie2
module load kallisto
module load bowtie
module load python3
module load java
cd /proj/sens2018579/TCR/tcr/tracer-master''', file=outfh)
			print('./tracer_python3 assemble -s Hsap %s %s /proj/sens2018579/TCR/tracer-out > %s 2> %s'%(' '.join(list(sorted(fastqs))), sample, os.path.join(o.folder_out, sample+'_tracer_stdout.txt'), os.path.join(o.folder_out, sample+'_tracer_stderr.txt')), file=outfh)
