from __future__ import division
import argparse, os, subprocess, gzip, sys, shutil, contextlib
from ConfigParser import SafeConfigParser
import taskmanager_lessthreaded as taskmanager

def find_fastq(path):
	for name in os.listdir(path):
		path2 = os.path.join(path, name)
		if os.path.isdir(path2):
			for fastq in find_fastq(path2): yield fastq
		elif any(name.endswith(ending) for ending in ('fastq.gz','fastq', 'fq.gz', 'fq')):
			yield path2

def staralign(conf, o, sample, mate1, mate2):
	samplefolder = os.path.join(o.processeddatafolder, sample)
	mkdir(samplefolder)
	outfolder = os.path.join(samplefolder, 'rnastar')
	mkdir(outfolder)
	if not outfolder.endswith('/'): outfolder += '/'
	addcmd = []
	if mate1.endswith('.gz'):
		addcmd += ['--readFilesCommand', 'zcat']
	if mate2:
		cmd = [conf.get('rnastar', 'path'), '--runThreadN', conf.get('rnastar', 'processes'), '--genomeDir', conf.get('rnastar', 'index'), '--readFilesIn', mate1, mate2]+addcmd+['--outFileNamePrefix', outfolder]
	else:
		cmd = [conf.get('rnastar', 'path'), '--runThreadN', conf.get('rnastar', 'processes'), '--genomeDir', conf.get('rnastar', 'index'), '--readFilesIn', mate1]+addcmd+['--outFileNamePrefix', outfolder]
	try: subprocess.check_call(cmd)
	except:
		raise Exception, 'Failed: '+' '.join(cmd)
	os.remove(os.path.join(outfolder, 'SJ.out.tab'))  # a file created by STAR

def staralign_viruses(conf, o, sample, mate1, mate2):
	samplefolder = os.path.join(o.processeddatafolder, sample)
	mkdir(samplefolder)
	outfolder = os.path.join(samplefolder, conf.get('viruses_rnastar', 'foldername'))
	mkdir(outfolder)
	if not outfolder.endswith('/'): outfolder += '/'
	addcmd = []
	if mate1.endswith('.gz'):
		addcmd += ['--readFilesCommand', 'zcat']
	if mate2:
		cmd = [conf.get('viruses_rnastar', 'path'), '--runThreadN', conf.get('viruses_rnastar', 'processes'), '--genomeDir', conf.get('viruses_rnastar', 'index'), '--readFilesIn', mate1, mate2]+addcmd+['--outFileNamePrefix', outfolder]
	else:
		cmd = [conf.get('viruses_rnastar', 'path'), '--runThreadN', conf.get('viruses_rnastar', 'processes'), '--genomeDir', conf.get('viruses_rnastar', 'index'), '--readFilesIn', mate1]+addcmd+['--outFileNamePrefix', outfolder]
	try: subprocess.check_call(cmd)
	except:
		raise Exception, 'Failed: '+' '.join(cmd)
	os.remove(os.path.join(outfolder, 'SJ.out.tab'))  # a file created by STAR

def ShannonH(L):
	import math
	from collections import defaultdict
	posD = defaultdict(int)
	for l in L:
		posD[l] += 1
	total_coverage = float(sum(posD.values()))
	coverage_fractions = [c/total_coverage for c in posD.values()]
	return -sum(c*math.log(c) for c in coverage_fractions)

def postprocess_mapcont_fromsorted(inputfile, outputfile, rpkmfile=None, sample=None, keepmore_red=True):
	# inputfile should be a coordinate-sorted sam file, including header
	import os
	from collections import defaultdict
	# perhaps do like https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0922-z and calculate Shannon entropy H for transcript coverage evenness, but recalulated as TIN=100*(e^H)/len, would benefit memory-consumption-wise from pre-sorting by position
	lengths = {}
	evenness = defaultdict(float)
	speciescounts = defaultdict(int)
	uniquepos = defaultdict(int)
	total = 0
	mapped = 0
	lastrseq = '*'
	speciesorder = []
	
	normreads = None
	totalreads = None
	if rpkmfile is not None:
		sample_i = 1
		with open(rpkmfile, 'rU') as infh:
			for line in infh:
				p = line.rstrip('\r\n').split('\t')
				if line.startswith('#samples'):
					if sample is not None:
						sample_i = p.index(sample)
				elif line.startswith('#normalizationreads'):
					normreads = float(p[sample_i])
				elif line.startswith('#allmappedreads'):
					totalreads = float(p[sample_i])
				elif not line.startswith('#'): break
	
	def add_item(locations, lastrseq):
		speciescounts[lastrseq] = len(locations)
		evenness[lastrseq] = 100*2.71828**ShannonH(locations)/lengths[lastrseq]
		uniquepos[lastrseq] = len(set(locations))
	
	with open(inputfile, 'r') as infh:
		for line in infh:
			if line.startswith('@SQ'):
				p = line.rstrip('\r\n').split()
				rseq = p[1].split(':')[1]
				lengths[rseq] = int(p[2].split(':')[1])
				speciesorder.append(rseq)
			if line.startswith('@'): continue
			p = line.rstrip('\r\n').split('\t')
			total += 1
			rseq = p[2]
			if rseq == '*': continue
			mapped += 1
			if rseq != lastrseq:
				if lastrseq != '*':
					add_item(locations, lastrseq)
				locations = []
				lastrseq = rseq
			locations.append(int(p[3]))
		if lastrseq != '*':
			add_item(locations, lastrseq)
	
	with open(outputfile, 'w') as outfh:
		outfh.write('\t'.join([str(totalreads if totalreads is not None else total), 'total reads'])+'\n')
		outfh.write('\t'.join([str(mapped), 'mapped reads'])+'\n')
		fields = ['reads', 'mapped_pos','TIN_score', 'ref_length', 'reference']
		if normreads: fields.append('RPKM')
		outfh.write('\t'.join(fields)+'\n')
		for species in speciesorder:
			fields = [str(speciescounts[species]), str(uniquepos[species]), str(evenness[species]), str(lengths[species]), species]
			if normreads: fields.append(str(speciescounts[species]*1e9/normreads/lengths[species]))
			outfh.write('\t'.join(fields)+'\n')
	try: os.chmod(outputfile, 0o774)
	except: pass
	
	# remove the sam file
	if not keepmore_red:
		try: os.remove(inputfile)
		except: pass

def run_rpkmforgenes(conf, o, sample, alignmentfile, outputfile):
	cmd = ['python2', conf.get('rpkmforgenes', 'path'), '-i', alignmentfile, '-n', sample, '-a', conf.get('rpkmforgenes', 'geneannotation'), '-o', outputfile] + conf.get('rpkmforgenes', 'options').split()
	try: subprocess.check_call(cmd)
	except:
		raise Exception, 'Failed: '+' '.join(cmd)

def run_cufflinks(conf, o, sample, alignmentfile, cufflinks_dir):
	cmd = [conf.get('cufflinks', 'path'), '-G', conf.get('cufflinks', 'annotationgtf'), '-o', cufflinks_dir, alignmentfile, '-p', conf.get('cufflinks', 'processes')] + conf.get('cufflinks', 'options').split()
	try: subprocess.check_call(cmd, stderr=open('/dev/null'))
	except:
		print cmd
		raise Exception, 'Failed: '+' '.join(cmd)

def run_fastqc(conf, o, sample, fastqfiles):
	samplefolder = os.path.join(o.processeddatafolder, sample)
	mkdir(samplefolder)
	outfolder = os.path.join(samplefolder, 'fastqc')
	mkdir(outfolder)
	if not outfolder.endswith('/'): outfolder += '/'
	cmd = [conf.get('fastqc', 'path')] + fastqfiles + ['--outdir='+outfolder, '--quite']
	#cmd = [conf.get('fastqc', 'path')] + fastqfiles + ['--outdir='+outfolder]
	subprocess.check_call(cmd)

def zcat_to_tmp(fastqs, cmd='zcat'):
	if not fastqs: return '', None
	import tempfile
	if fastqs.split(',')[0].endswith('.fastq.gz') and cmd == 'cat':
		suffix = '.fastq.gz'
	elif fastqs.split(',')[0].endswith('.fastq') and cmd == 'zcat':
		suffix = '.fastq'
		cmd = 'cat'
	else: suffix = '.fastq'
	tmpfq = tempfile.mkstemp(suffix=suffix)[1]
	with open(tmpfq, 'w') as outfh:
		pr = subprocess.Popen([cmd]+fastqs.split(','), stdout=outfh)
	return tmpfq, pr

def run_mitcr(conf, o, sample, fastqfiles, outfiles):
	if not fastqfiles: raise Exception
	mate1,pr1 = zcat_to_tmp(','.join(fastqfiles), cmd='cat')
	if pr1.wait(): raise Exception
	try:
		samplefolder = os.path.join(o.processeddatafolder, sample)
		mkdir(samplefolder)
		cmd1 = ['java', '-Xmx2g', '-jar', conf.get('mitcr', 'path')]+ conf.get('mitcr', 'options').split() +['-gene', 'TRA', mate1, outfiles[0]]
		cmd2 = ['java', '-Xmx2g', '-jar', conf.get('mitcr', 'path')]+ conf.get('mitcr', 'options').split() +['-gene', 'TRB', mate1, outfiles[1]]
		if o.debugmode:
			print ' '.join(cmd1)
			print ' '.join(cmd2)
			subprocess.check_call(cmd1)
			subprocess.check_call(cmd2)
		else:
			with open('/dev/null', 'w') as nullfh:
				subprocess.check_call(cmd1, stderr=nullfh)
				subprocess.check_call(cmd2, stderr=nullfh)
	finally:
		if mate1: os.remove(mate1)

def run_mixcr(conf, o, sample, mate1_fqs, mate2_fqs, outfile):
	mixcr_path = conf.get('mixcr', 'path')
	samplefolder = os.path.join(o.processeddatafolder, sample)
	mkdir(samplefolder)
	outfolder = os.path.join(samplefolder, 'mixcr')
	mkdir(outfolder)
	if len(mate1_fqs) != 1 or len(mate2_fqs) != 1: raise Exception
	cmd = range(4)
	cmd[0] = [mixcr_path, 'align', '-p', 'rna-seq', '-r', os.path.join(outfolder, 'log1.txt'), mate1_fqs[0], mate2_fqs[0], os.path.join(outfolder, 'alignments.vdjca')]
	cmd[1] = [mixcr_path, 'assemblePartial', os.path.join(outfolder, 'alignments.vdjca'), os.path.join(outfolder, 'alignment_contigs.vdjca')]
	cmd[2] = [mixcr_path, 'assemble', '-r', os.path.join(outfolder, 'log2.txt'), os.path.join(outfolder, 'alignment_contigs.vdjca'), os.path.join(outfolder, 'clones.clns')]
	cmd[3] = [mixcr_path, 'exportClones', os.path.join(outfolder, 'clones.clns'), outfile]
	for c in cmd:
		if o.debugmode:
			print ' '.join(c)
		subprocess.check_call(c)

def star_bam(conf, o, sample, samfile, bamfile):
	assert bamfile.endswith('.bam')
	if not os.path.exists(samfile): raise Exception
	cmd1 = ['samtools', 'view', '-Sb', samfile] # -q 255 filters for unique alignments
	pr1=subprocess.Popen(cmd1, stdout=subprocess.PIPE)
	cmd2 = ['samtools', 'sort', '-']
	with open(bamfile, 'wb') as outfh:
		subprocess.check_call(cmd2, stdin=pr1.stdout, stdout=outfh)
	cmd3 = ['samtools', 'index', bamfile]
	subprocess.check_call(cmd3)
	if pr1.wait(): raise Exception

def sortsam(samfile, conf):
	if not os.path.exists(samfile):
		print samfile
		raise Exception
	tmpfile = samfile+'_sorted.sam'
	cmd = ['java', '-jar', conf.get('viruses_rnastar', 'picard_sortsam_path'), 'I='+samfile, 'O='+tmpfile, 'SORT_ORDER=coordinate', 'QUIET=true', 'VERBOSITY=WARNING']
	subprocess.check_call(cmd)
	shutil.move(tmpfile, samfile)

def bowtiealign_hla(conf, o, sample, mate1, mate2, zipped=False):
	if not mate1: raise Exception
	if zipped:
		mate1,pr1 = zcat_to_tmp(mate1)
		mate2,pr2 = zcat_to_tmp(mate2)
		if pr1.wait(): raise Exception
		if mate2 and pr2.wait(): raise Exception
	try:
		samplefolder = os.path.join(o.processeddatafolder, sample)
		mkdir(samplefolder)
		outfile = os.path.join(samplefolder, sample+'_bowtiealignHLA.sam')
		if mate2:
			cmd = [conf.get('bowtie_hla', 'path'), '-p', conf.get('bowtie_hla', 'processes'), '-n', '3', '--sam', conf.get('bowtie_hla', 'index'), '-1', mate1, '-2', mate2, outfile]
		else:
			cmd = [conf.get('bowtie_hla', 'path'), '-p', conf.get('bowtie_hla', 'processes'), '-n', '3', '--sam', conf.get('bowtie_hla', 'index'), mate1, outfile]
		subprocess.check_call(cmd)
	finally:
		if zipped:
			if mate1: os.remove(mate1)
			if mate2: os.remove(mate2)

def htseq_count(file_in, file_out, annotation_file, options_str, htseqcount_path):
	with open(file_out, 'w') as outfh:
		subprocess.check_call([htseqcount_path, '-f', 'bam', '-r', 'pos']+options_str.split()+[file_in, annotation_file], stdout=outfh)

def bowtie2(conf, o, sample, mate1, mate2, zipped, samfile, confsection='bowtie2'):
	if not mate1: raise Exception
	if zipped:
		mate1,pr1 = zcat_to_tmp(mate1)
		mate2,pr2 = zcat_to_tmp(mate2)
		if pr1.wait(): raise Exception
		if mate2 and pr2.wait(): raise Exception
	try:
		samplefolder = os.path.join(o.processeddatafolder, sample)
		mkdir(samplefolder)
		outfile = samfile
		if mate2:
			cmd = [conf.get('bowtie', 'path'), '-p', conf.get('bowtie', 'processes'), '-n', '3', '--sam', conf.get('bowtie', 'index'), '-1', mate1, '-2', mate2, outfile]
		else:
			cmd = [conf.get('bowtie', 'path'), '-p', conf.get('bowtie', 'processes'), '-n', '3', '--sam', conf.get('bowtie', 'index'), mate1, outfile]
		subprocess.check_call(cmd)
	finally:
		if zipped:
			if mate1: os.remove(mate1)
			if mate2: os.remove(mate2)

def get_read_length(fastqfile):
	openfunc = gzip.open if fastqfile.endswith('.gz') else open
	infh = openfunc(fastqfile, 'r')
	try:
		for li, line in enumerate(infh):
			if li == 1:
				return len(line.rstrip())
	finally:
		infh.close()

def run_seq2hla(conf, o, sample, mate1, mate2, zipped, seq2hla_file):
	if not mate1: raise Exception
	readlength = get_read_length(mate1.split(',')[0])
	if zipped:
		mate1,pr1 = zcat_to_tmp(mate1)
		mate2,pr2 = zcat_to_tmp(mate2)
		if pr1.wait(): raise Exception
		if mate2 and pr2.wait(): raise Exception
	current_wd = os.getcwd()
	try:
		samplefolder = os.path.join(o.processeddatafolder, sample)
		mkdir(samplefolder)
		outfolder = os.path.abspath(os.path.join(samplefolder, 'seq2hla'))
		mkdir(outfolder)
		if mate2:
			cmd = ['python2', os.path.join(conf.get('seq2hla', 'folderpath'), 'seq2HLA.py'), '-1', mate1, '-2', mate2, '-l', str(readlength), '-r', os.path.join(outfolder, sample)]
		else:
			raise Exception
		os.chdir(conf.get('seq2hla', 'folderpath'))
		with open(seq2hla_file, 'w') as outfh:
			subprocess.check_call(cmd, stdout=outfh)
	finally:
		if zipped:
			if mate1: os.remove(mate1)
			if mate2: os.remove(mate2)
		os.chdir(current_wd)

def get_normreads(rpkmfile, sample, allmapped=False):
	with open(rpkmfile, 'rU') as infh:
		for line in infh:
			if not line.startswith('#'): break
			p = line.rstrip().split('\t')
			if p[0] == '#samples': samples = p[1:]
			elif p[0] == ('#allmappedreads' if allmapped else '#normalizationreads'): normreads = p[1:]
	i = samples.index(sample)
	return normreads[i]

def rpkms_hla(conf, o, sample, samfile):
	samplefolder = os.path.join(o.processeddatafolder, sample)
	mkdir(samplefolder)
	outfile = os.path.join(samplefolder, sample+conf.get('rpkms_hla', 'suffix'))
	rpkmfile = os.path.join(o.processeddatafolder, sample, sample+conf.get('rpkmforgenes', 'suffix'))
	allmapped = conf.getboolean('rpkms_hla', 'allmappednorm')
	normreads = get_normreads(os.path.join(o.processeddatafolder, sample, sample+conf.get('rpkmforgenes', 'suffix')), sample, allmapped)
	cmd = ['python2', conf.get('rpkms_hla', 'path'), samfile, '-D', normreads] + conf.get('rpkms_hla', 'options').split()
	with open(outfile, 'w') as outfh:
		subprocess.check_call(cmd, stdout=outfh)

def bowtie_and_rpkm_HLA(conf, o, sample, mate1, mate2, zipped, samfile):
	if (not o.dontredo) or (not os.path.exists(samfile)):
		bowtiealign_hla(conf, o, sample, mate1, mate2, zipped)
	rpkms_hla(conf, o, sample, samfile)
	if os.path.exists(samfile):
		os.remove(samfile)

def mkdir(folder):
	if not os.path.exists(folder):
		os.mkdir(folder)

def list_samples(o):
	if o.samples is None:
		if o.align or o.fastqc or o.mitcr or o.viruses or o.mixcr:
			return os.listdir(o.rawdatafolder)
		else:
			return os.listdir(o.processeddatafolder)
	else:
		return o.samples

def star_rmsam(conf, o, sample, alignmentfile, bamfile):
	samexists = os.path.exists(alignmentfile)
	bamexists = os.path.exists(bamfile)
	sortingbam = os.path.exists(bamfile[:-4]+'.0000.bam')
	if sortingbam: raise Exception, 'Encountered temporary sorting bam file for %s'%sample
	if samexists and bamexists and not sortingbam:
		os.remove(alignmentfile)

def get_paired_fastq(sample, o):
	fastqfiles = list(sorted(find_fastq(os.path.join(o.rawdatafolder,sample))))
	mate1 = ','.join([f for f in fastqfiles if f.endswith('_1.fastq.gz') or f.endswith('_R1_001.fastq') or f.endswith('_R1_001.fastq.gz')])
	mate2 = ','.join([f for f in fastqfiles if f.endswith('_2.fastq.gz') or f.endswith('_R2_001.fastq') or f.endswith('_R2_001.fastq.gz')])
	if not mate2: mate1 = ','.join(fastqfiles)
	if all(not fq.endswith('.gz') for fq in fastqfiles):
		zipped = False
	elif all(fq.endswith('.gz') for fq in fastqfiles):
		zipped = True
	else:
		raise Exception
	return mate1, mate2, zipped

def find_alignmentfile(sample, o):
	alignmentfile = os.path.join(o.processeddatafolder, sample, 'rnastar', sample+'.bam')
	if not os.path.exists(alignmentfile) and not o.makebam:
		alignmentfile = os.path.join(o.processeddatafolder, sample, 'rnastar', 'Aligned.out.sam')
	return alignmentfile

@contextlib.contextmanager
def cd(newdir):
#https://stackoverflow.com/questions/431684/how-do-i-cd-in-python/24176022#24176022
    prevdir = os.getcwd()
    os.chdir(newdir)
    try:
        yield
    finally:
        os.chdir(prevdir)

def vdjpuzzle(sample, o, conf): # untested
	outfolder = os.path.join(o.processeddatafolder, sample, 'vdjpuzzle')
	mkdir(outfolder)
	with cd(outfolder):
		with open(os.path.join(outfolder, 'processlog.txt')) as logfh:
			subprocess.check_call(['bash', conf.get('VDJPuzzle', 'path'), os.path.join(o.rawdatafolder, sample)], stdout=logfh, stderr=logfh) # shouldn't work, wrong folder structure

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-r', '--rawdatafolder', help='e.g. J.Lundeberg_14_22', required=True)
	parser.add_argument('-a', '--processeddatafolder', required=True)
	parser.add_argument('-F', '--fastqc', action='store_true')
	parser.add_argument('-A', '--align', action='store_true')
	parser.add_argument('-B', '--makebam', action='store_true', help=argparse.SUPPRESS)
	parser.add_argument('--rmsamifbam', action='store_true', help=argparse.SUPPRESS)
	parser.add_argument('-R', '--rpkms', action='store_true')
	parser.add_argument('--HTSeqcount', action='store_true')
	parser.add_argument('--cufflinks', action='store_true')
	parser.add_argument('-T', '--mitcr', action='store_true', help=argparse.SUPPRESS)
	parser.add_argument('-X', '--mixcr', action='store_true')
	parser.add_argument('-H', '--hla', action='store_true', help=argparse.SUPPRESS)
	parser.add_argument('-S', '--seq2hla', action='store_true')
	parser.add_argument('-V', '--viruses', action='store_true')
	parser.add_argument('--alignbowtie2', action='store_true')
	parser.add_argument('-s', '--samples', nargs='+')
	parser.add_argument('--debugmode', action='store_true')
	parser.add_argument('--dontredo', action='store_true')
	parser.add_argument('--conf')
	parser.add_argument('--processes', type=int, help=argparse.SUPPRESS) #redundant
	o = parser.parse_args()
	if o.align:
		o.makebam = True
		o.rmsamifbam = True
	
	if o.conf is None:
		o.conf = os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])), 'configuration_process_rnaseq.ini')
	
	conf = SafeConfigParser()
	if not conf.read(o.conf): raise Exception
	
	tasklist = taskmanager.Tasklist(o.processes if o.processes is not None else int(conf.get('general', 'processes')), verbose=True, singleprocess=o.debugmode)
	
	if o.fastqc:
		mkdir(o.processeddatafolder)
		for sample in list_samples(o):
			fastqfiles = list(sorted(find_fastq(os.path.join(o.rawdatafolder,sample))))
			if o.dontredo and os.path.exists(os.path.join(o.processeddatafolder, sample, 'fastqc')): continue
			tasklist.add(run_fastqc, (conf, o, sample, fastqfiles), sample=sample, group='fastqc', num_p=1, timeout_hours=1)
	
	if o.align:
		mkdir(o.processeddatafolder)
		for sample in list_samples(o):
			fastqfiles = list(sorted(find_fastq(os.path.join(o.rawdatafolder,sample))))
			mate1 = ','.join([f for f in fastqfiles if f.endswith('_1.fastq.gz')])
			mate2 = ','.join([f for f in fastqfiles if f.endswith('_2.fastq.gz')])
			alignmentfile = os.path.join(o.processeddatafolder, sample, 'rnastar', 'Aligned.out.sam')
			bamfile = os.path.join(o.processeddatafolder, sample, 'rnastar', sample+'.bam')
			if not mate2: mate1 = ','.join(fastqfiles)
			if o.dontredo and (os.path.exists(alignmentfile) or os.path.exists(bamfile)): continue
			print mate1, mate2
			tasklist.add(staralign, (conf, o, sample, mate1, mate2), sample=sample, group='staralign', num_p=int(conf.get('rnastar', 'processes')), timeout_hours=16, maxingroup=1)
	
	if o.makebam:
		for sample in list_samples(o):
			alignmentfile = os.path.join(o.processeddatafolder, sample, 'rnastar', 'Aligned.out.sam')
			bamfile = os.path.join(o.processeddatafolder, sample, 'rnastar', sample+'.bam')
			if o.dontredo and os.path.exists(bamfile): continue
			tasklist.add(star_bam, (conf, o, sample, alignmentfile, bamfile), sample=sample, group='star_bam', num_p=1, timeout_hours=1, waitfor=['staralign'])
	
	if o.rmsamifbam:
		for sample in list_samples(o):
			alignmentfile = os.path.join(o.processeddatafolder, sample, 'rnastar', 'Aligned.out.sam')
			bamfile = os.path.join(o.processeddatafolder, sample, 'rnastar', sample+'.bam')
			tasklist.add(star_rmsam, (conf, o, sample, alignmentfile, bamfile), sample=sample, group='star_rmsam', num_p=1, timeout_hours=1, waitfor=['staralign', 'star_bam'])
	
	if o.rpkms:
		for sample in list_samples(o):
			alignmentfile = find_alignmentfile(sample, o)
			rpkmfile = os.path.join(o.processeddatafolder, sample, sample+conf.get('rpkmforgenes', 'suffix'))
			if o.dontredo and os.path.exists(rpkmfile): continue
			tasklist.add(run_rpkmforgenes, (conf, o, sample, alignmentfile, rpkmfile), sample=sample, group='rpkms', num_p=1, timeout_hours=10, waitfor=['staralign', 'star_bam'], maxingroup=20)
	
	if o.cufflinks:
		for sample in list_samples(o):
			alignmentfile = find_alignmentfile(sample, o)
			cufflinks_dir = os.path.join(o.processeddatafolder, sample, 'cufflinks')
			mkdir(cufflinks_dir)
			if o.dontredo and os.path.exists(os.path.join(cufflinks_dir, 'isoforms.fpkm_tracking')): continue
			tasklist.add(run_cufflinks, (conf, o, sample, alignmentfile, cufflinks_dir), sample=sample, group='cufflinks', num_p=int(conf.get('cufflinks', 'processes')), timeout_hours=2, waitfor=['staralign', 'star_bam'], maxingroup=20)
	
	if o.HTSeqcount:
		for sample in list_samples(o):
			alignmentfile = find_alignmentfile(sample, o)
			countfile = os.path.join(o.processeddatafolder, sample, sample+conf.get('htseqcount', 'suffix'))
			if o.dontredo and os.path.exists(countfile): continue
			tasklist.add(htseq_count, (alignmentfile, countfile, conf.get('htseqcount', 'annotationgtf'), conf.get('htseqcount', 'options'), conf.get('htseqcount', 'path')), sample=sample, group='htseqcount', num_p=1, timeout_hours=1, waitfor=['staralign', 'star_bam'])
	
	if o.viruses:
		mkdir(o.processeddatafolder)
		for sample in list_samples(o):
			fastqfiles = list(sorted(find_fastq(os.path.join(o.rawdatafolder,sample))))
			foldername = conf.get('viruses_rnastar', 'foldername')
			mate1 = ','.join([f for f in fastqfiles if f.endswith('_1.fastq.gz')])
			mate2 = ','.join([f for f in fastqfiles if f.endswith('_2.fastq.gz')])
			alignmentfile = os.path.join(o.processeddatafolder, sample, foldername, 'Aligned.out.sam')
			print alignmentfile
			bamfile = os.path.join(o.processeddatafolder, sample, foldername, sample+'.bam')
			if not mate2: mate1 = ','.join(fastqfiles)
			if o.dontredo and (os.path.exists(alignmentfile) or os.path.exists(bamfile)): pass
			else:
				tasklist.add(staralign_viruses, (conf, o, sample, mate1, mate2), sample=sample, group='viruses_rnastar', num_p=int(conf.get('viruses_rnastar', 'processes')), timeout_hours=16, maxingroup=1)
			counts_file = os.path.join(o.processeddatafolder, sample, foldername, 'counts.txt')
			if o.dontredo and os.path.exists(counts_file): pass
			else:
				tasklist.add(sortsam, (alignmentfile,conf), sample=sample, group='viruses_sortsam', num_p=int(conf.get('viruses_rnastar', 'processes')))
				if conf.getboolean('viruses_rnastar', 'include_RPKM'):
					rpkmfile = os.path.join(o.processeddatafolder, sample, sample+conf.get('rpkmforgenes', 'suffix'))
					tasklist.add(postprocess_mapcont_fromsorted, (alignmentfile, counts_file, rpkmfile, sample), waitfor=['viruses_rnastar', 'viruses_sortsam', 'rpkms'], group='viruses_postprocess', sample=sample)
				else:
					tasklist.add(postprocess_mapcont_fromsorted, (alignmentfile, counts_file), waitfor=['viruses_rnastar', 'viruses_sortsam', 'rpkms'], group='viruses_postprocess', sample=sample)
	
	if o.mitcr:
		mkdir(o.processeddatafolder)
		for sample in list_samples(o):
			fastqfiles = list(sorted(find_fastq(os.path.join(o.rawdatafolder,sample))))
			mitcrfiles = [os.path.join(o.processeddatafolder, sample, sample+'_TRA'+conf.get('mitcr', 'suffix')), os.path.join(o.processeddatafolder, sample, sample+'_TRB'+conf.get('mitcr', 'suffix'))]
			if o.dontredo and all(os.path.exists(mitcrfile) for mitcrfile in mitcrfiles): continue
			tasklist.add(run_mitcr, (conf, o, sample, fastqfiles, mitcrfiles), sample=sample, group='mitcr', num_p=4, timeout_hours=10, waitfor=['rpkms'])
	
	if o.mixcr:
		if o.debugmode: print 'mixcr'
		mkdir(o.processeddatafolder)
		for sample in list_samples(o):
			fastqfiles = list(sorted(find_fastq(os.path.join(o.rawdatafolder,sample))))
			mate1 = [f for f in fastqfiles if f.endswith('_1.fastq.gz') or '_R1_' in os.path.split(f)[-1]]
			mate2 = [f for f in fastqfiles if f.endswith('_2.fastq.gz') or '_R2_' in os.path.split(f)[-1]]
			mixcrfile = os.path.join(o.processeddatafolder, sample, 'mixcr', 'clones.txt')
			if o.debugmode: print sample, mate1, mate2, fastqfiles
			if o.dontredo and os.path.exists(mixcrfile): continue
			tasklist.add(run_mixcr, (conf, o, sample, mate1, mate2, mixcrfile), sample=sample, group='mixcr', num_p=1, timeout_hours=30, waitfor=[])
	
	if o.hla:
		mkdir(o.processeddatafolder)
		for sample in list_samples(o):
			mate1, mate2, zipped = get_paired_fastq(sample, o)
			samfile = os.path.join(o.processeddatafolder, sample, sample+'_bowtiealignHLA.sam')
			rpkmfileHLA = os.path.join(o.processeddatafolder, sample, sample+conf.get('rpkms_hla', 'suffix'))
			if o.dontredo:
				if os.path.exists(rpkmfileHLA): continue
			tasklist.add(bowtie_and_rpkm_HLA, (conf, o, sample, mate1, mate2, zipped, samfile), sample=sample, group='bowtiealign_hla', num_p=int(conf.get('bowtie_hla', 'processes')), waitfor=['rpkms'])
	
	if o.alignbowtie2:
		mkdir(o.processeddatafolder)
		for sample in list_samples(o):
			mate1, mate2, zipped = get_paired_fastq(sample, o)
			samfile_b = os.path.join(o.processeddatafolder, sample, sample+conf.get('bowtie2', 'suffix'))
			if o.dontredo:
				if os.path.exists(samfile_b): continue
			tasklist.add(bowtie2, (conf, o, sample, mate1, mate2, zipped, samfile_b, 'bowtie2'), sample=sample, group='bowtie2align', num_p=int(conf.get('bowtie2', 'processes')))
	
	if o.seq2hla:
		mkdir(o.processeddatafolder)
		for sample in list_samples(o):
			mate1, mate2, zipped = get_paired_fastq(sample, o)
			seq2hla_file = os.path.abspath(os.path.join(o.processeddatafolder, sample, 'seq2hla', 'progress.txt'))

			if o.dontredo and os.path.exists(seq2hla_file): continue
			tasklist.add(run_seq2hla, (conf, o, sample, mate1, mate2, zipped, seq2hla_file), sample=sample, group='seq2hla', num_p=6, timeout_hours=3)
	
	tasklist.waitforall()
