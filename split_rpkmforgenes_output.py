import argparse, sys, time
from collections import defaultdict

def loadlist(infile, index=None, func=None, ignore='#', ignorelines=0):
	""" returns array of strings """
	infileh = open(infile, "rU")
	for i in range(ignorelines):
		infileh.readline()
	outarray = [l.strip().split()[0] for l in infileh.readlines()]
	infileh.close()
	if ignore is not None:
		outarray = [l for l in outarray if not l.startswith(ignore)]
	if index is not None:
		outarray = [l.split("\t")[index] for l in outarray]	
	if func is not None:
		outarray = list(map(func, outarray))
	return outarray

def name_match(sample, group, o):
	if o.name_match == 'part': return any(part in sample for part in group)
	elif o.name_match == 'exact': return sample in group

if '__main__' == __name__:
	opts = argparse.ArgumentParser()
	opts.add_argument('inf', help='from rpkmforgenes.py')
	opts.add_argument('outf', nargs='?', default='/dev/stdout')
	opts.add_argument('-i', '--include', nargs='+', metavar='sample')
	opts.add_argument('-e', '--exclude', nargs='*', default=[], metavar='sample')
	opts.add_argument('-g', '--genelist', metavar='genelistfile')
	opts.add_argument('-G', '--genes', nargs='+', metavar='gene')
	opts.add_argument('-T', '--tableformat', action='store_true')
	opts.add_argument('-T2', '--tableformat2', action='store_true')
	opts.add_argument('-T3', '--tableformat3', action='store_true')
	opts.add_argument('-R', '--tableformatreads', action='store_true')
	opts.add_argument('-C', '--csvformat', action='store_true')
	opts.add_argument('-i2', '--include_really', nargs='+', default=[], help=argparse.SUPPRESS)
	opts.add_argument('-s', '--samplelistfile')
	opts.add_argument('--ignoremissingsamples', action='store_true')
	opts.add_argument('-r', '--rename', nargs=2, action='append', default=[], metavar='name')
	opts.add_argument('-m', '--name_match', choices=['part', 'exact'], default='part')
	o = opts.parse_args()
	
	if o.tableformatreads:
		if o.tableformat: raise Exception
		o.tableformat = True
	if o.tableformat2:
		#if o.tableformat: raise Exception
		o.tableformat = True
	if o.tableformat3:
		#if o.tableformat or o.tableformat2: raise Exception
		o.tableformat = True
	
	genelist = None
	if o.genelist is not None:
		genelist = set(loadlist(o.genelist))
	if o.genes is not None:
		if genelist is None: genelist = set()
		genelist |= set(o.genes)
	
	namecounter = defaultdict(int)
	
	with open(o.inf, 'rU') as infh:
		with open(o.outf, 'w') as outfh:
			for line in infh:
				p = line.rstrip('\r\n').split('\t')
				if p[0] == '#samples':
					sample_i = [i for i,s in enumerate(p[1:]) if ((o.include is None or name_match(s, o.include, o)) and not name_match(s, o.exclude, o) or name_match(s, o.include_really, o))]
					num_samples = len(p[1:])
					
					if o.samplelistfile:
						# those both in sample_i and in the list file
						# if no o.include and no o.exclude, sample_i will have all samples
						with open(o.samplelistfile, 'rU') as samplefileh:
							samples_in_list = [line.rstrip('\r\n') for line in samplefileh]
						sample_i_set = set(sample_i)
						if o.ignoremissingsamples:
							sample_i_in_list = [p[1:].index(s) for s in samples_in_list if s in p[1:]]
						else:
							sample_i_in_list = [p[1:].index(s) for s in samples_in_list] # raises ValueError if rpkmfile is missing any sample
						sample_i = [i for i in sample_i_in_list if i in sample_i_set]
					
					for pair in o.rename:
						for si, sample in enumerate(p):
							if si == 0: continue
							if sample == pair[0]:
								p[si] = pair[1]
					
					if o.tableformat:
						if o.tableformat3:
							print >>outfh, 'gene\tIDs\t'+'\t'.join('%s'%p[i+1] for i in sample_i)
						elif o.tableformat2:
							print >>outfh, 'gene\t'+'\t'.join('%s'%p[i+1] for i in sample_i)
						else:
							print >>outfh, 'gene\t'+'\t'.join('"%s"'%p[i+1] for i in sample_i)
					elif o.csvformat:
						print >>outfh, 'sample,'+','.join(p[i+1] for i in sample_i)
					else:
						try:
							print >>outfh, '\t'.join(p[i+1] for i in [-1]+sample_i)
						except:
							print p, sample_i, len(p)
							raise
				elif p[0].startswith('#arguments'):
					if not o.tableformat and not o.csvformat:
						if p[0] == '#arguments':
							print >>outfh, '\t'.join(p+[' '.join(sys.argv), 'time: ' + time.asctime()])
						else:
							print >>outfh, '\t'.join(['#arguments']+[p[0][10:]]+p[1:]+[' '.join(sys.argv), 'time: ' + time.asctime()])
				elif p[0] in ('#allmappedreads', '#genemappedreads', '#normalizationreads'):
					if not o.tableformat and not o.csvformat:
						try:
							print >>outfh, '\t'.join(p[i+1] for i in [-1]+sample_i)
						except IndexError:
							while len(p) <= max(sample_i) +1:
								p.append('0')
							print >>outfh, '\t'.join(p[i+1] for i in [-1]+sample_i)
				else:
					if genelist is not None and \
					 not any(gene in genelist for gene in p[0].split('+')) and \
					 not any(gene in genelist for gene in p[1].split('+')):
						continue
					if o.tableformat:
						name1 = p[0]
						name2 = p[1]
						if name2 in ('.','+','-',''):
							name = name1
						elif name1 in ('.','+','-',''):
							name = name2
						elif o.tableformat3:
							name = name1.replace('+','|') + '\t' + name2.replace('+','|')
						elif o.tableformat2:
							name = name1.replace('+',' ')
						else:
							name = name1.replace('+',' ')+' '+name2
						if not o.tableformatreads:
							print >>outfh, name + '\t' + '\t'.join(p[i+2] for i in sample_i)
						else:
							print >>outfh, name + '\t' + '\t'.join(p[i+num_samples+2] for i in sample_i)
					elif o.csvformat:
						name = p[0]
						namecounter[name] += 1
						if namecounter[name] > 1:
							name = '%s(%d)'%(name, namecounter[name])
						print >>outfh, name + ',' + ','.join(p[i+2] for i in sample_i)
					else:
						if len(p) - 2 == 2*num_samples:
							print >>outfh, '\t'.join([p[i+2] for i in [-2, -1]+sample_i] + [p[i+num_samples+2] for i in sample_i])
						else:
							print >>outfh, '\t'.join(p[i+2] for i in [-2, -1]+sample_i)
