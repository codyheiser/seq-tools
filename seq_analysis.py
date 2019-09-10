# -*- coding: utf-8 -*-
"""
seq_analysis.py
@author: Cody Heiser


usage: seq_analysis.py [-h] [--truncate] sequence

Analyze a DNA sequence

positional arguments:
  sequence    input .txt file or string with DNA sequence

optional arguments:
  -h, --help  show this help message and exit
  --truncate  truncate printing of complimentary sequences to reduce console output

"""

import os
import argparse
import zipfile
import gzip
import pandas as pd

class dna_sequence:

	def __init__(self, seq):
		self.sequence = seq # string containing DNA sequence
		self.counts = {
		'A' : self.sequence.count('A')+self.sequence.count('a'),
		'T' : self.sequence.count('T')+self.sequence.count('t'),
		'C' : self.sequence.count('C')+self.sequence.count('c'),
		'G' : self.sequence.count('G')+self.sequence.count('g')
		}
		self.seq_length = sum(self.counts.values()) # total length of sequence in bp


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
        'A' : {'DNA':'T', 'RNA':'U'},
		'T' : {'DNA':'A', 'RNA':'A'},
		'G' : {'DNA':'C', 'RNA':'C'},
		'C' : {'DNA':'G', 'RNA':'G'},
		'a' : {'DNA':'T', 'RNA':'U'},
		't' : {'DNA':'A', 'RNA':'A'},
		'g' : {'DNA':'C', 'RNA':'C'},
		'c' : {'DNA':'G', 'RNA':'G'},
        ' ' : {'DNA':' ', 'RNA':' '},
        '\t' : {'DNA':'\t', 'RNA':'\t'},
        '\n' : {'DNA':'\n', 'RNA':'\n'},
        ',' : {'DNA':',', 'RNA':','}
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
			return round((self.counts['A'] + self.counts['T'])*2 +
			(self.counts['G'] + self.counts['C'])*4, 3)

		else:
			return round(64.9 + 41*(self.counts['G'] + self.counts['C'] - 16.4)/(self.counts['A'] +
				self.counts['T'] +
				self.counts['G'] +
				self.counts['C']), 3)


	def mass_calc(self):
		'''
		Perform calculation of molar mass of oligo sequence, assuming synthesized, anhydrous ssDNA
		'''
		return round(self.counts['A']*313.21 + self.counts['T']*304.2 +
			self.counts['C']*289.18 + self.counts['G']*329.21 - 61.96, 3)


	@classmethod
	def from_file(cls, datafile):
		'''initialize object from outside file (datafile)'''
		filetype = os.path.splitext(datafile)[1] # extract file extension to save as metadata

		if filetype == '.zip': # if compressed, open the file and update filetype
			zf = zipfile.ZipFile(datafile)
			datafile = zf.open(os.path.splitext(datafile)[0]) # update datafile with zipfile object
			filetype = os.path.splitext(os.path.splitext(datafile)[0])[1] # update filetype

		elif filetype == '.gz': # if file is g-zipped, read accordingly
			f = gzip.open(datafile, 'r') # open file containing DNA sequence. 'r' = read only
			data = f.read() # read DNA sequence into variable
			f.close() # close file after reading

		else:
			f = open(datafile, 'r') # open file containing DNA sequence. 'r' = read only
			data = f.read() # read DNA sequence into variable
			f.close() # close file after reading

		return cls(data)



if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Analyze a DNA sequence')
	parser.add_argument('sequence', help='input .txt file or string with DNA sequence')
	parser.add_argument('--truncate', help='truncate printing of complimentary sequences to reduce console output',
		action='store_true')
	parser.add_argument('--prettyprint', help='print full sequence metrics to console output',
		action='store_true')
	args = parser.parse_args()

	# if 'sequence' argument is a file, read it into a string
	if os.path.isfile(args.sequence):
		s = dna_sequence.from_file(args.sequence)

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
		print_len = len(s.sequence)+1
		post_script = ''

	if args.prettyprint:
		# print stats about the sequence to the console
		print("\nSequence:\n    3`-{}{}-5`".format(s.sequence[:print_len], post_script))
		print("\nReverse compliment:\n    3`-{}{}-5`".format(s.rev_comp('DNA')[:print_len], post_script))
		print("\nRNA reverse compliment:\n    3`-{}{}-5`".format(s.rev_comp('RNA')[:print_len], post_script))
		print("\nLength: {} bp".format(s.seq_length))
		print("\nBase counts: {}".format(s.counts))
		print("\nGC content: {} %".format(s.gc_content()))
		print("\nMelting temperature (Tm): {} C".format(s.tm_calc()))
		print("\nMolar mass: {} g/mol\n\n".format(s.mass_calc()))

	else:
		print(s.rev_comp('DNA')[:print_len], post_script)
