#!/usr/bin/python

#==================================================
#  Created on: Aug 2017
#      Author: xiao
#==================================================

#import libraries
import os, sys
import subprocess
import time 

#==================================================
#
# Usage: 
#	'python run_MSA.py <latinname_seqfiles_dir_path> <MSA_program_exe_path>' 
#	OR
#	'nohup python run_MSA.py <latinname_seqfiles_dir_path> <MSA_program_exe_path>'
#
#==================================================

#Defined functions


def run_MSA(latinname_seqfiles_dir_path, MSA_program_exe_path_argv):
	#MSA_program_exe_path = '/home/xiao/MUSCLE/muscle3.8.31_i86linux64'
	MSA_program_exe_path = MSA_program_exe_path_argv
	#MSA_opt1 = " -seqtype protein "
	MSA_opt2 = " -in "
	MSA_opt3 = " -fastaout "
#	MSA_opt4 = ' -maxiters 1 -diags1 -sv -distance1 kbit20_3'

	for seq in os.listdir(latinname_seqfiles_dir_path):
		if '.fasta' in seq and 'latin' in seq and 'aligned_' not in seq:
			if os.access(latinname_seqfiles_dir_path, os.R_OK):
				os.chdir(latinname_seqfiles_dir_path)
				
				gene_name = seq.strip('.fasta').split('$')[1]
				input_name = seq.split('$')[0] + '\$' + seq.split('$')[1]

				output_name = 'aligned_' + gene_name + '.fasta'

				cmd = MSA_program_exe_path  + MSA_opt2 + input_name + MSA_opt3 + output_name
				print(cmd)
				#run MSA as input
				os.system(cmd)

				print("finish MSA " + gene_name + '\n')				
		
#==================================================

if __name__ == '__main__':
	start_time = time.time()
	sys.stdout = open('MSA_report.txt','wr')

	latinname_seqfiles_dir_path = sys.argv[1]
	MSA_program_exe_path_argv = sys.argv[2]
	
	run_MSA(latinname_seqfiles_dir_path, MSA_program_exe_path_argv)

	print("Time Used: --- %s seconds ---" % (time.time() - start_time) + '\n')
	sys.stdout.close()	
