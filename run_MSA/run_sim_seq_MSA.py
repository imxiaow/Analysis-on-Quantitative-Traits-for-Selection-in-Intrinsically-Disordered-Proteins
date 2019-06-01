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
#	'python run_sim_MSA.py <simulated_seqfiles_dir_path> <MSA_program_exe_path>' 
#	OR
#	'nohup python run_sim_MSA.py <simulated_seqfiles_dir_path> <MSA_program_exe_path>'
#
#==================================================

#Defined functions

def run_sim_MSA(sim_dir_path, MSA_program_exe_path_argv):
	#MSA_program_exe_path = '/home/xiao/MUSCLE/muscle3.8.31_i86linux64'
	MSA_program_exe_path = MSA_program_exe_path_argv
	MSA_opt2 = " -in "
	MSA_opt3 = " -fastaout "
	
	for sim_seq in os.listdir(sim_dir_path):
		if 'aa.fasta' in sim_seq and '_idr_' in sim_seq and 'sim_' in sim_seq and 'aligned' not in sim_seq:
			if os.access(sim_dir_path, os.R_OK):
				os.chdir(sim_dir_path)

				input_name = sim_seq
				output_name = 'aligned_' + sim_seq
			
				cmd = MSA_program_exe_path  + MSA_opt2 + input_name + MSA_opt3 + output_name
				print(cmd)
				#run MSA as input
				os.system(cmd)
				print("finish MSA " + sim_seq + '\n')

#==================================================

if __name__ == '__main__':
	start_time = time.time()
	sys.stdout = open('sim_MSA_report.txt','wr')

	path_name = sys.argv[1]
	MSA_program_exe_path_argv = sys.argv[2]

	run_sim_MSA(path_name,MSA_program_exe_path_argv)

	
	print("Time Used: --- %s seconds ---" % (time.time() - start_time) + '\n')
	sys.stdout.close()	