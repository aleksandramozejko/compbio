#!/usr/bin/python

from __future__ import division
from Bio import SeqIO
from os.path import join
from os import listdir
import sys
import argparse
import math
import random
import itertools

########## Functions' definitions ##########
def seq_parser(f):
	'''Given a fasta file f returns a generator of sequences from the file'''
	return SeqIO.parse(open(f), 'fasta')

def family_size(gen):
	'''Given a generator gen of sequences from a family returns an integer
	representing the size of the family.
	'''
 	return len(list(gen))

def sample(family, size):
	'''Given a family of protein sequences and its size returns a list of 
	of sequences sampled from the family. If the size of a family is greater
	than 20 the number of sequences selected from the family equals
	a ceiling of 0.1 times the size of the family. Otherwise it returns all 
	sequences from the family.
	'''
	vect = [x for x in xrange(size)]
	if size > 10:
		sample_size = int(math.ceil(0.1 * size))
	else:
		sample_size = size
	sample_indices = sorted(random.sample(vect, sample_size))
	sample = []
	counter = 0
	for protein in family:
	 	if counter in sample_indices:
	 		sample.append(str(protein.seq))
		counter += 1
	return sample

def create_pairs(lst):
	'''Given a list lst returns a list pairs of all elements of the list
	lst paired (without repetitions).
	'''
	pairs = []
	for pair in itertools.combinations(lst, 2):
		pairs.append(list(pair))
	return pairs

def parse_matrix(matrix_file):
	'''Reads substitution matrix from matrix_file and returns that matrix
	represented as a two-dimensional list matrix.
	'''
	with open(matrix_file) as f:
		matrix = [line.strip().split() for line in f if '#' not in line]
	return matrix

def score(l1, l2, subs_matrix):
	'''Given two strings l1 and l2 representing nucleid acids and 
	a substitution matrix returns an int representing a score of l1 to l2
	amino acid transformation, which is determined by the substitution matrix.
	'''
	nucleic_acids = subs_matrix[0]
	l1_ind = nucleic_acids.index(l1)
	l2_ind = nucleic_acids.index(l2)
	return int(subs_matrix[l1_ind + 1][l2_ind + 1])

def find_max(matrix):
	'''Given a two-dimensional list matrix returns its maximal element.
	'''
	return max([max(x) for x in matrix])

def alignment_score(seqa, seqb, subs_matrix):
	'''Given two strings seqa, seqb (representing sequences) and 
	a two-dimensional list subs_matrix returns the score of an local
	alignment with no gaps of seqa and seqb. 
	'''
	score_matrix = [[0] * (len(seqb) + 1) for i in xrange(len(seqa) + 1)]
	for i in xrange(1,len(seqa)+1):
		for j in xrange(1,len(seqb)+1):
			single_score = score(seqa[i-1], seqb[j-1], subs_matrix)
			score_matrix[i][j] = max(0, score_matrix[i-1][j-1] + single_score)
	return find_max(score_matrix)

def family_score(sample, subs_matrix):
	'''Given a list sample from the family and a substitution matrix
	subs_matrix represented as a list of lists returns the score
	of the alignment for the whole family, i.e. arithmetic mean
	of single scores of all pairwise alignments from sample.
	'''
	pairs = create_pairs(sample)
	scores = []
	for pair in pairs:
		scores.append(alignment_score(pair[0], pair[1], subs_matrix))
	return sum(scores)/len(scores)

########## Parsing arguments from the command line ##########
parser = argparse.ArgumentParser()
parser.add_argument('-d', 
					metavar = 'DIR', required = True,
 					help ='Input directories with PAM/BLOSUM matrices')
parser.add_argument('-f',
					metavar = 'DIR',
					help ='Input directory with fasta files')
parser.add_argument('--filenames',
					nargs = '+',
					help = 'Input fasta files')
parser.add_argument('--results',
					metavar = 'FILE',
					help = 'Output file')
args = parser.parse_args()


matdir = args.d
matdir_lst = listdir(matdir)
results = open(args.results, 'w')

# If directory with the fasta files is given
if args.f:
	famdir = args.f
	famdir_lst = listdir(famdir)
	for fam in famdir_lst:
		# Printing the results to stdout
		results.write('family name: ' + fam)
		family = seq_parser(join(famdir, fam))
		fam_size = family_size(seq_parser(join(famdir, fam)))
		results.write('family size: ' + str(fam_size) + '\n')
		samp = sample(family, fam_size)
		for mat in matdir_lst:
			matrix = parse_matrix(join(matdir, mat))
	 		results.write(' '.join([matdir, mat, 'score:', str(family_score(samp, matrix))]))
	 	results.write('\n')

#If single fasta files are given
if args.filenames:
	filenames = args.filenames
	for fil in filenames:
		#Printing the results to stdout
		results.write('family name: ' + fil)
		family = seq_parser(fil)
		fam_size = family_size(seq_parser(fil))
		results.write('family size: ' + str(fam_size) + '\n')
		samp = sample(family, fam_size)
		for mat in matdir_lst:
			matrix = parse_matrix(join(matdir, mat))
	 		results.write(' '.join([matdir, mat, 'score:', str(family_score(samp, matrix))]))
	 	results.write('\n')




