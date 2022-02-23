from __future__ import division
import argparse, os, subprocess, sys
from ConfigParser import SafeConfigParser
import taskmanager_lessthreaded as taskmanager

def find_fastq(path):
	for name in os.listdir(path):
		path2 = os.path.join(path, name)
		if os.path.isdir(path2):
			for fastq in find_fastq(path2): yield fastq
		elif any(name.endswith(ending) for ending in ('fastq.gz','fastq')):
			yield path2

def list_samples(o):
	if o.samples is None:
		return os.listdir(o.folder_in)
	else:
		return o.samples

def mkdir(folder):
	if not os.path.exists(folder):
		os.mkdir(folder)

def trim_galore(fastqs_in, samplefolder_out, conf):
	try: options = conf.get('trim_galore', 'options').split()
	except: options = []
	mkdir(samplefolder_out)
	subprocess.check_call([conf.get('trim_galore', 'path')] + fastqs_in + options + ['-o', samplefolder_out], stderr=open('/dev/null','w'))

# command e.g. python ~/src/run_trim_galore.py --debugmode INBOX/P8654/ trimmed_fastq/ -s P8654_1012 --samplenametranslation samplenametranslation.txt

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('folder_in')
	parser.add_argument('folder_out')
	parser.add_argument('-s' ,'--samples', nargs='+')
	parser.add_argument('--conf')
	parser.add_argument('--debugmode', action='store_true')
	parser.add_argument('--dontredo', action='store_true')
	parser.add_argument('--samplenametranslation', metavar='tabdelfile')
	parser.add_argument('--processes', type=int, help=argparse.SUPPRESS)
	o = parser.parse_args()
	
	if o.conf is None:
		o.conf = os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])), 'configuration_process_rnaseq.ini')
	conf = SafeConfigParser()
	if not conf.read(o.conf): raise Exception
	
	tasklist = taskmanager.Tasklist(o.processes if o.processes is not None else int(conf.get('general', 'processes')), verbose=True, singleprocess=o.debugmode)
	
	samplenametranslation = dict()
	if o.samplenametranslation:
		with open(o.samplenametranslation, 'rU') as infh:
			for line in infh:
				p = line.rstrip('\r\n').split('\t')
				samplenametranslation[p[0]] = p[1]
	
	mkdir(o.folder_out)
	for sample in list_samples(o):
		samplefolder_in = os.path.join(o.folder_in, sample)
		fastqs_in = list(find_fastq(samplefolder_in))
		samplefolder_out = os.path.join(o.folder_out, samplenametranslation.get(sample, sample))
		if o.dontredo and os.path.exists(samplefolder_out): continue
		tasklist.add(trim_galore, (fastqs_in, samplefolder_out, conf), sample=sample, group='trim_galore')
	
	tasklist.waitforall()
