#!/usr/bin/python

#==================================================
#  Created on: Aug 2017
#      Author: xiao
#==================================================

#import libraries
import sys
import os
import time

#==================================================
#
# Usage: 
#	'python change_names.py <file_dirpath>'
#
#==================================================

def read_write_file(data_file, file_name):
	new_filename = 'human_' + file_name.split('$')[1]
	write_filename = new_filename
	write_file = open(write_filename, 'wr')
	for line in data_file:
		write_file.write(line)
	write_file.close()

#==================================================

if __name__ == '__main__':
	path_name = sys.argv[1]

	for filename in os.listdir(path_name):
		if 'human' in filename and '.fasta' in filename:
			os.chdir(path_name)
			f = open(filename,'r')
			read_write_file(f, filename)