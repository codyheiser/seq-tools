# -*- coding: utf-8 -*-
"""
seq_analysis.py
@author: Cody Heiser


usage: Heiser_7.py [-h] [--truncate] sequence

Analyze a DNA sequence

positional arguments:
  sequence    input .txt file or string with DNA sequence

optional arguments:
  -h, --help  show this help message and exit
  --truncate  truncate printing of complimentary sequences to reduce console output

"""

import os
import argparse

class dna_sequence:

	def __init__(self, seq):
		self.sequence = seq.upper() # string containing DNA sequence
		self.seq_length = len(seq) # length of sequence in bp
		self.counts = {
		'A':self.sequence.count('A'),'T':self.sequence.count('T'),
		'C':self.sequence.count('C'),'G':self.sequence.count('G')
		}
		# if letters other than A, T, C, and G are given, raise error
		if self.seq_length > self.counts['A']+self.counts['T']+self.counts['G']+self.counts['C']:
			raise ValueError('Sequence must contain only nucleotide base letters (A, T, C, G)')

	def prettyprint_metrics(self):
		'''
		print base metrics to console.
		'''
		print('\nSequence has the following base counts:\n\tA : {}\n\tT : {}\n\tG : {}\n\tC : {}\n'.format(self.counts['A'], self.counts['T'], self.counts['G'], self.counts['C']))

	def rev_comp(self, nucleic_acid):
		'''
		Generate reverse compliment to given sequence, s,
		in desired nucleic acid format (i.e. DNA or RNA).
		Output is reported 3' to 5'.
		'''
		# define base complements
		bases = {
		'A':{'DNA':'T', 'RNA':'U'},
		'T':{'DNA':'A', 'RNA':'A'},
		'G':{'DNA':'C', 'RNA':'C'},
		'C':{'DNA':'G', 'RNA':'G'}
		}
		# return reverse complimentary string (3' - 5')
		return ''.join([bases[x][nucleic_acid] for x in self.sequence])[::-1]

	def gc_content(self):
		'''
		Calculate percent GC content of sequence, s
		'''
		return round(((self.counts['G'] + self.counts['C'])/self.seq_length)*100, 3)

	def tm_calc(self):
		'''
		Perform basic calculation of melting temperature (Tm) based on nucleotide sequence
		'''
		if self.seq_length < 14:
			return round((self.sequence.count('A') + self.sequence.count('T'))*2 +
			(self.sequence.count('G') + self.sequence.count('C'))*4, 3)

		else:
			return round(64.9 + 41*(self.sequence.count('G') + self.sequence.count('C') - 16.4)/(self.sequence.count('A') +
				self.sequence.count('T') +
				self.sequence.count('G') +
				self.sequence.count('C')), 3)

	def mass_calc(self):
		'''
		Perform calculation of molar mass of oligo sequence, assuming synthesized, anhydrous ssDNA
		'''
		return round(self.sequence.count('A')*313.21 + self.sequence.count('T')*304.2 +
			self.sequence.count('C')*289.18 + self.sequence.count('G')*329.21 - 61.96, 3)



if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Analyze a DNA sequence')
	parser.add_argument('sequence', help='input .txt file or string with DNA sequence')
	parser.add_argument('--truncate', help='truncate printing of complimentary sequences to reduce console output',
		action='store_true')
	args = parser.parse_args()

	# if 'sequence' argument is a file, read it into a string
	if os.path.isfile(args.sequence):
		# get sequence data to play with
		f = open(args.sequence, 'r') # open file containing DNA sequence. 'r' = read only
		seq = f.read() # read DNA sequence into variable
		f.close() # close file after reading

	# otherwise, just use provided string
	else:
		seq = args.sequence

	# create dna_sequence object with provided sequence
	s = dna_sequence(seq)

	# if user provided flag to truncate output for large sequence, provide print limit of 80 characters
	if args.truncate and s.seq_length > 80:
		print_len = 81
		post_script = '...'

	# otherwise, return entire strings for reverse complements
	else:
		print_len = s.seq_length+1
		post_script = ''

	# print stats about the sequence to the console
	print("\nSequence:\n    3`-{}{}-5`".format(s.sequence[:print_len], post_script))
	print("\nReverse compliment:\n    3`-{}{}-5`".format(s.rev_comp('DNA')[:print_len], post_script))
	print("\nRNA reverse compliment:\n    3`-{}{}-5`".format(s.rev_comp('RNA')[:print_len], post_script))
	print("\nLength: {} bp".format(s.seq_length))
	print("\nBase counts: {}".format(s.counts))
	print("\nGC content: {} %".format(s.gc_content()))
	print("\nMelting temperature (Tm): {} C".format(s.tm_calc()))
	print("\nMolar mass: {} g/mol\n\n".format(s.mass_calc()))
