#! /usr/bin/python3

import argparse, os, sys
from Bio import SeqIO

def makelen2(array, empty):
	if len(array) > 2: raise Exception
	while(len(array)) < 2: array.append(empty)

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('tracer_fasta', nargs='+')
	o = parser.parse_args()
	
	print('\t'.join(['sample', 'TRA_1_seq', 'TRA_2_seq', 'TRB_1_seq', 'TRB_2_seq', 'TRA_number', 'TRB_number', 'TRA_1_descr', 'TRA_1_V', 'TRA_1_J', 'TRA_2_descr', 'TRA_2_V', 'TRA_2_J', 'TRB_1_descr', 'TRB_1_V', 'TRB_1_J', 'TRB_2_descr', 'TRB_2_V', 'TRB_2_J']))
	for fastafile in o.tracer_fasta:
		sample = os.path.split(fastafile)[-1].rsplit('.f',1)[0].rsplit('_TCR',1)[0]
		
		with open(fastafile, 'r') as infh:
			A = []; Adescr = []
			B = []; Bdescr = []
			for seqitem in SeqIO.parse(infh, 'fasta'):
				name, sequence = seqitem.id, str(seqitem.seq)
				descr = name.split('|')[-1]
				descr = [descr, descr.rsplit('_',2)[0], descr.split('_')[-1]]
				if '|TCR|A|' in name:
					A.append(sequence)
					Adescr.append(descr)
				elif '|TCR|B|' in name:
					B.append(sequence)
					Bdescr.append(descr)
				else: raise Exception
		An = len(A); Bn = len(B)
		makelen2(A, ''); makelen2(B, ''); makelen2(Adescr, ['','','']); makelen2(Bdescr, ['','',''])
		print('\t'.join([sample]+A+B+[str(An), str(Bn)]+Adescr[0]+Adescr[1]+Bdescr[0]+Bdescr[1]))
