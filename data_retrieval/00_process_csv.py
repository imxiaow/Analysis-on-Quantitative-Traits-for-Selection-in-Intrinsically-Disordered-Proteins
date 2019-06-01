#!/usr/bin/python

#==================================================
#  Created on: Aug 2017
#      Author: xiao
#==================================================

#import libraries
import sys
import os
import re, csv

#==================================================
#
# This script is use to process csv file that download from uniprot web, which it must contain column 'Gene names (primary)'
# This script is reading csv files and extracting primary gene names, writing into a txt file.
# It will also checking cases where multiple symbols are assigned for same protein, and record this information in multisymbolforsameprotein file.
# Place this script at the same directory of csv file.
# Enter csv filename by user input, for example: uniprot-all.csv
#
#==================================================

#Defined function

def read_data(data_file, data_filename):
	reader = csv.DictReader(data_file)
	write_filename = 'multisymbolforsameprotein_' + data_filename.split('.csv')[0] +'.txt'
	symbols = []
	multisymbol = []

	for row in reader:
		gene_symbols = row['Gene names  (primary )'].split(' ')
		if (len(gene_symbols)>1):
			multisymbol.append(gene_symbols)
			for symbl_ in gene_symbols:
				if ';' in symbl_:
					symbl_ = symbl_.strip(';')
				symbols.append(symbl_)

		if (len(gene_symbols)==1):
			for sym in gene_symbols:
				if ';' in sym:
					sym = sym.strip(';')
				symbols.append(sym)

		else:
			pass

	if len(multisymbol) != 0:
		write_file = open(write_filename,'wr')
		write_file.write(str(multisymbol) +'\n')
		write_file.close()

	return symbols

#==================================================

if __name__ == '__main__':
	csv_filename = raw_input('Please enter the name of csv file: \n')
	write_filename = csv_filename.split('.csv')[0] + '_genes.txt'
	write_file = open(write_filename,'wr')
	csv_file= open(csv_filename ,'r')
	print('reading your csv file..')
	symbols = read_data(csv_file, csv_filename)

	for symbol in symbols:
		if symbol != '':
			write_file.write(symbol+'\n')
	write_file.close()
	print('finish :)')
