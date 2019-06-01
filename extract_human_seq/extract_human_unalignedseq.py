#!/usr/bin/python

#==================================================
#  Created on: Aug 2017
#      Author: xiao
#==================================================

#import libraries
import os
import sys
import time
#==================================================
#
# Usage: 
#	'python extract_human_unalignedseq.py <file_dirpath_need_to_extract>'
#
#==================================================

#Defined functions
def find_human_seq_wf(gene_name, path):
	with open(path) as data_file:
		write_file_name = 'human_unalignedseq$' + gene_name +'.fasta'
		write_file = open(write_file_name,'wr')
		flag = False
		for line in data_file:
			if line.startswith('>')and 'homo_sapiens' in line:
				flag = True
			if line.startswith('>') and 'homo_sapiens' not in line:
				flag = False
			if flag is True:
				write_file.write(line)
		write_file.close()
		print("Finish extract human sequence of " + gene_name + '\n')
	data_file.close()


#==================================================

if __name__ == '__main__':
	sys.stdout = open('extract_humanseq_report.txt','wr')
	pathname = sys.argv[1]

	for subsubfiles in os.listdir(pathname):
		if '.fasta' in subsubfiles and 'latin_names' in subsubfiles:
			if os.access(os.path.join(pathname, subsubfiles), os.R_OK):
				os.chdir(pathname)
				gene_name = subsubfiles.strip('.fasta').split('$')[1]
				print("extract human protein from " + gene_name + ' now\n')				
				find_human_seq_wf(gene_name,(os.path.join(pathname, subsubfiles)))

	sys.stdout.close()
