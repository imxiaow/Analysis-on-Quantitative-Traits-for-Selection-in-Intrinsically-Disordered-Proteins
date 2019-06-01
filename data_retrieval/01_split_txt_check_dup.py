#!/usr/bin/python

#==================================================
#  Created on: Aug 2017
#      Author: xiao
#==================================================

#import libraries
import sys
import os
from itertools import izip_longest
import shutil

#==================================================
#
#  This program will ask about user preferences such as number of genes in each splitted files.
#  Users need to enter the path to target file, for example: /Users/.../work_dir/target_file.txt
#  And then enter choice "y" or "n" for changing default numbers or not.
#  After all the user input, this program will split genes into separate file and create new directories corresponding to file name,
#  and then move matched files into its corresponded directory, this prepare you ready for data retrieval program. (seperated dirtories)
#  Also, it will check duplicated gene symbols in each splitted file in each directory, and recorded,
#  this is convenient for you to check total number of genes later after you run data retrieval program.
#
#==================================================

#Defined functions and variables 

n = 2000 #default setting of number in each splitted file

def grouper(n, iterable, fillvalue=None): #helper function use to split gene txt file
	'''
	collect data into fixed length chucks or blocks
	group(3, 'ABCDEFG', 'x') ---> ABC DEF GXX
	'''
	args = [iter(iterable)] * n 

	return izip_longest(fillvalue= fillvalue, *args)


def check_duplitcate_gene_name(gene_file_name, dir_path): #function use to check duplicated gene name in a splitted file
	gene_file = open(dir_path + gene_file_name,'r')

	gene_lst = []
	for line in gene_file:
		line = line.split('\n')[0]
		gene_lst.append(line)

	checked = []
	duplicated = []
	for genes in gene_lst:
		if genes not in checked:
			checked.append(genes)
		else:
			duplicated.append(genes)

	if len(duplicated)!= 0: # there are duplicated gene names in this file
		print('duplicated gene symbol found, recording..')
		duplicated_file = open(dir_path + 'duplicatedgene_'+ gene_file_name, 'w')
		duplicated_file.write(gene_file_name + '\n' + 'duplicated the following genes in gene_file\n'+str(duplicated))
		duplicated_file.close()

	else: 
		print('no duplicated gene symbol :)')


#==================================================

if __name__ == '__main__':
	file_wantto_split_path = raw_input('Please enter the path of text_file that needed to run split: \n')
	#get user input
	while not (os.access(file_wantto_split_path, os.R_OK) and '.txt' in file_wantto_split_path): #user input is not valid
		file_wantto_split_path = raw_input('Please re-enter the path of text_file that needed to run split: ')

	else: #user input is valid
		print('Found file path :)')
		check_number = raw_input('The default number in each splited file is 2000.\nIf you dont want to change it, enter: n\nIf you would like to change it, enter: y\nYour choice(y/n): \n')
		#checking user input choice
		while not (check_number == 'y' or check_number =='n'): #user input is not valid
			check_number= raw_input('Please enter your choice(y/n): \n')

		else: # valid user input
			if check_number == 'y': #user input number
				number = raw_input('Please enter the number of genes you perfer in each splited files: \n')
				number = int(number)

				with open(file_wantto_split_path) as f:
					for i, g in enumerate(grouper(number, f, fillvalue= ''), 1): #split gene names in file into separated files
						with open('subfile_{0}.txt'.format(i * number), 'w')as fout: 
							fout.writelines(g) 

						current_dir = os.getcwd()
						new_dir = current_dir + '/' + ('subfile_{0}.txt'.format(i * number)).split('.txt')[0] + '/'
						if not os.path.exists(new_dir): #making new directory that correspond to splited gene files
							os.makedirs(new_dir)

				for txt_file in os.listdir(os.getcwd() + '/'):
					if 'subfile_' in txt_file and '.txt' in txt_file:
						source = os.getcwd() + '/' + txt_file
						destination = os.getcwd() + '/' + txt_file.split('.txt')[0] + '/' + txt_file
						shutil.move(source, destination) #moving corresponding gene file to its own named directory
						source2 = os.getcwd() + '/' + 'main.py'
						destination2 = os.getcwd() + '/' + txt_file.split('.txt')[0] + '/'+ 'main.py'
						shutil.copy2(source2, destination2)
						source3 = os.getcwd() + '/' + 'fetch_fetch_fetch.py'
						destination3 = os.getcwd() + '/' + txt_file.split('.txt')[0] + '/'+ 'fetch_fetch_fetch.py'
						shutil.copy2(source3, destination3)
						source4 = os.getcwd() + '/' + 'check_geneRefID_by_UniprotID.py'
						destination4 = os.getcwd() + '/' + txt_file.split('.txt')[0] + '/'+ 'check_geneRefID_by_UniprotID.py'
						shutil.copy2(source4, destination4)
						source5 = os.getcwd() + '/' + 'extract_high_conf_geneRefIDs.py'
						destination5 = os.getcwd() + '/' + txt_file.split('.txt')[0] + '/''extract_high_conf_geneRefIDs.py'
						shutil.copy2(source5, destination5)

				for item in os.listdir(os.getcwd() + '/'): 
					path_n = os.getcwd() + '/' + item
					if os.path.isdir(path_n):
						for gene_file in os.listdir(path_n + '/'):
							if 'subfile_' in gene_file and '.txt' in gene_file and 'duplicated' not in gene_file:
								if os.access(path_n+'/'+ gene_file, os.R_OK):
									print('checking duplicated gene_name in ' + gene_file + '..')
									check_duplitcate_gene_name(gene_file, (path_n+'/')) #checking duplicated gene names in each subfile directory

				print('finished :)')

			else: #default n = 2000
				with open(file_wantto_split_path) as f:
					for i, g in enumerate(grouper(n, f, fillvalue= ''), 1):
						with open('subfile_{0}.txt'.format(i * n), 'w')as fout:
							fout.writelines(g)

						current_dir = os.getcwd()
						new_dir = current_dir + '/' + ('subfile_{0}.txt'.format(i * n)).split('.txt')[0] + '/'
						if not os.path.exists(new_dir):
							os.makedirs(new_dir)

				for txt_file in os.listdir(os.getcwd() + '/'):
					if 'subfile_' in txt_file and '.txt' in txt_file:
						source = os.getcwd() + '/' + txt_file
						destination = os.getcwd() + '/' + txt_file.split('.txt')[0] + '/' + txt_file
						shutil.move(source, destination)
						source2 = os.getcwd() + '/' + 'main.py'
						destination2 = os.getcwd() + '/'+ txt_file.split('.txt')[0] + '/' + 'main.py'
						shutil.copy2(source2, destination2)
						source3 = os.getcwd() + '/' + 'fetch_fetch_fetch.py'
						destination3 = os.getcwd() + '/' + txt_file.split('.txt')[0] + '/'+ 'fetch_fetch_fetch.py'
						shutil.copy2(source3, destination3)
						source4 = os.getcwd() + '/' + 'check_geneRefID_by_UniprotID.py'
						destination4 = os.getcwd() + '/' + txt_file.split('.txt')[0] + '/'+ 'check_geneRefID_by_UniprotID.py'
						shutil.copy2(source4, destination4)
						source5 = os.getcwd() + '/' + 'extract_high_conf_geneRefIDs.py'
						destination5 = os.getcwd() + '/' + txt_file.split('.txt')[0] + '/'+ 'extract_high_conf_geneRefIDs.py'
						shutil.copy2(source5, destination5)


				for item in os.listdir(os.getcwd() + '/'):
					path_n = os.getcwd() + '/' + item
					if os.path.isdir(path_n):
						for gene_file in os.listdir(path_n + '/'):
							if 'subfile_' in gene_file and '.txt' in gene_file and 'duplicated' not in gene_file:
								if os.access(path_n+'/'+ gene_file, os.R_OK):
									print('checking duplicated gene_name in ' + gene_file + '..')
									check_duplitcate_gene_name(gene_file, (path_n+'/'))

				print('finished :)')