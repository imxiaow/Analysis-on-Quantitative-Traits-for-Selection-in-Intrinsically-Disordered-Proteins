#!/usr/bin/python

import os
import sys
import time
import subprocess

if __name__ == '__main__':
	pathname = sys.argv[1]
	#DISO_program_exe_path = '/home/xiao/DISOPRED/run_disopred.pl '
	DISO_program_exe_path = sys.argv[2]

	for subfile in os.listdir(pathname):
		if 'human' in subfile and '.fasta' in subfile:

			cmd = DISO_program_exe_path + ' ' + pathname + subfile
			print(cmd)

			os.system(cmd)
			print("finish Disordered prediction " + subfile)
