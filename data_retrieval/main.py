#!/usr/bin/python

#======================================================================
#  Created on: Aug 2017
#      Author: xiao
#======================================================================

#import libraries
import requests, sys
import json
import pprint
import xml.etree.ElementTree as ET 
import check_geneRefID_by_UniprotID as CKGENEID #helper functions
import fetch_fetch_fetch as FFFROOT #helper functions
import extract_high_conf_geneRefIDs as EXTRACT #helper functions
import os
import time
import shutil

#======================================================================
#
# This script uses user prereference, based on commandline argument, it will check for duplication if not yet checked before.
# Standard outputs are redirected to report file. 
# High confidence orthologues sequence along with human reference sequence are written in all_seqs_with_latin_names$gene_symbol.fasta,
# gene name wih less than 11 orthologues does not have fasta file, gene name is recorded as txt file.
# protein sequence IDs are also recorded in its corresponded file, all_seqs_with_protIDs$gene_symbol.fasta,
# duplicated species in each gene_names fasta files are checked and recorded as duplicated_species_genename.txt
# failed cases are already handled in helper functions that from imported scripts.
#
# The speed of data retrieving is heavily depend on network speed and data size, also, there is rate limiting on request times.
# Rate limiting is handled by waiting, however, it will reach maximum depth resursion if waited too many times.
# Be aware about the number of processes you are using, especially on static router server.
#
# Usage: 
#	'python main.py <dir_path> <gene_txtfile> <y/n>' 
#	OR 
#	'nohup python main.py <dir_path> <gene_txtfile> <y/n>'
#
#======================================================================

#Defined variables 


server = "https://rest.ensembl.org" #REST Server
find_orthologues = "/homology/symbol/human/{}?" #REST API
symbol_look_up="/xrefs/name/human/{}" #REST API
look_up_iDs = "/xrefs/id/{}?" #REST API
get_sequence = "/sequence/id/{}?type=protein" #REST API

#Defined functions

def extract_Names_ProtID(symbol, root, geneRefId_Lst_H): #get dictionary of species name with its protein id
	names_ProtId_dict = {}
	geneRef_Id_homo = {}
	prot_Id_homo = []

	for c in root.findall('{http://orthoXML.org/2011/}species'): #parsing xml files
		for c1 in c.findall('{http://orthoXML.org/2011/}database'):
			for c2 in c1.findall('{http://orthoXML.org/2011/}genes'):
				for gene in c2:
					if gene.attrib['id'] in geneRefId_Lst_H:
						value_lst = []
						if c.attrib['name'] not in names_ProtId_dict:
							names_ProtId_dict[c.attrib['name']] = value_lst
							names_ProtId_dict[c.attrib['name']].append(gene.attrib['protId'])

						else:
							names_ProtId_dict[c.attrib['name']].append(gene.attrib['protId'])

						if (c.attrib['name'] == 'homo_sapiens'):
							if gene.attrib['geneId'] not in geneRef_Id_homo:
								geneRef_Id_homo[gene.attrib['geneId']] = gene.attrib['id']
								prot_Id_homo.append(gene.attrib['protId'])

	print(geneRef_Id_homo)

	if len(geneRef_Id_homo) < 1: 
		print("no human reference ID")

	if len(geneRef_Id_homo) == 1: #exactly one gene reference ID for human
		pass

	if len(geneRef_Id_homo) > 1: #more than one gene reference ID for human, means synonym gene names are included, need to distinguish
		print("There are multiple human geneRefIDs!!(synonyms might included)\nChecking out real symbol GeneRefID...\n")
		geneRef_Id_homo_lst = []
		for key__ in geneRef_Id_homo:
			geneRef_Id_homo_lst.append(key__)

		check_info = CKGENEID.get_real_symbol_geneID(symbol, geneRef_Id_homo_lst,server, symbol_look_up, look_up_iDs)
		print(check_info)

		#checking the results from get_real_symbol_geneID, all failed cases already handled in check_geneRefID_by_UniprotID.py
		if ((check_info == -2) or (check_info == 1234) or (check_info == 5678) or (check_info == 9000) or (check_info == 91234) or (check_info == 95678)):
			return 11111

		else: # distinguished, exactly one matched gene_reference ID 
			acutal_query_human_geneRefID = check_info[0]
			print("re-extract names_ProtId_dict based on acutal_query_human_geneId...\n")
			for id__ in geneRef_Id_homo:
				if id__ == acutal_query_human_geneRefID:
					human_geneID = geneRef_Id_homo[id__]
			#extracting high confidence gene reference ID, human ref ID included
			high_conf_geneReflist_WithHuman = EXTRACT.re_extract_high_conf_geneRefID_H(human_geneID,root)
			#based on gene_ref ID, extrating latin name and protein ID
			names_ProtId_dict = EXTRACT.re_extract_Names_ProtID(root, high_conf_geneReflist_WithHuman)

	#--------------------------------------------------
	#checking duplicated species in high confidence file, recorded 
	count = 0
	duplicated_species_name = []
	for name in names_ProtId_dict:
		if(len(names_ProtId_dict[name])) > 1:
			duplicated_species_name.append(name)

		count += len(names_ProtId_dict[name])

	if len(duplicated_species_name)>0: #duplicated species found
		print("Attention!! Duplicated species(more than one sequence from same species): " + str(duplicated_species_name) + '\n')
		write_duplicatefilename = 'duplicated_speices_'+symbol + '.txt'
		write_duplicatefile = open(write_duplicatefilename, 'wr')
		write_duplicatefile.write(symbol +'\n'+'Duplicated species:\n'+str(duplicated_species_name))
		write_duplicatefile.close()
	#--------------------------------------------------

	print("The length of high_conf protID dictionary(human & duplicated species included) is : "+ str(count) +"\n")

	return names_ProtId_dict


def check_ortho_file(root):	#checking xml ortholog file before parsing
	check = 0
	c = root.find('{http://orthoXML.org/2011/}groups')

	if c is None: #when ortholog file is empty
		print("Empty orthologue file, may be a psudogenes!!!\n")
		return -1

	c1 = c.findall('{http://orthoXML.org/2011/}property')

	if c1 is None: #when ortholog file is empty
		print("Empty orthologue file, may be a psudogenes!!!\n")
		return -1

	homo_count = 0 #checking for number of human 
	for s in root.findall('{http://orthoXML.org/2011/}species'):
		if (s.attrib['name'] == 'homo_sapiens'):
			homo_count += 1

	if homo_count == 0: #when there is no human information in ortho file, can't retrieve anything.
		return -1

	return check 


def main_operations(symbol,path):
	command = symbol
	print("Welcome :) \n\nfinding orthologues from human reference protein of "+ command +" ...")
	root = FFFROOT.write_process_file_getroot(command,server)

	if root == None: #no information from ortholog file, already handled in other script.
		pass

	else:
		check = check_ortho_file(root)

		if check == -1: #check if the ortholog file is empty, if the file is empty, it may be a psudogenes(no proteins); or there's no human information 
			write_filename = 'psudogene_' + symbol + '.txt'
			write_file = open(write_filename,'wr')
			write_file.write(symbol + '\n' + "empty orthologue file, check if it's psudogene on Ensembl web." +'\n')
			write_file.close()

		if check == 0: #ortholog file is not empty
			print("finding high confidence orthologues from all orthologues ...\n")
			high_conf_geneReflist_WithHuman = EXTRACT.get_high_conf_geneRefID(EXTRACT.extract_high_conf_paired_geneRefID(root))
			
			if (len(high_conf_geneReflist_WithHuman) < 1): #recording information in case there is no high confidence orthologues at all
				write_filename = 'Failed_nohighor_' + symbol + '.txt'
				write_file = open(write_filename,'wr')
				write_file.write("No get_high_conf_geneRefID avaliable in this gene file. Failed\n")
				write_file.close()

			elif (len(high_conf_geneReflist_WithHuman) >= 1) and (len(high_conf_geneReflist_WithHuman) < 12): #recording information in case there is less than 12 orthologues sequences available
				write_filename = 'notenough_' + symbol + '.txt'
				write_file = open(write_filename,'wr')
				write_file.write(symbol + "Number of high_conf_gene is less than 12\n")
				write_file.close()

			else: #contain more than 12 orthologues 
				high_conf_LatinName_ProtID_Dic_WithHuman = extract_Names_ProtID(command, root, high_conf_geneReflist_WithHuman)
				if (high_conf_LatinName_ProtID_Dic_WithHuman != 11111) and (high_conf_LatinName_ProtID_Dic_WithHuman is not None): #success case
					print("getting sequences from high confidence orthologues ProtID ...\n")
					FFFROOT.write_sequence_file(high_conf_LatinName_ProtID_Dic_WithHuman,command) #extract protein sequence based on high confidence dictionary

					for file_name in os.listdir(path):#find seq_file with prot ID
						if (('.fasta' in file_name) and ('all_seqs_with_protIDs' in file_name ) and (command in file_name)):
							file_name_c = file_name.strip('.fasta').split('$')[1]

							if file_name_c == command:
								seq_file = open(file_name,'r')#open it for read
								#converting protein ID to corresponding latin name
								FFFROOT.process_seq_file_to_latinName_wt(high_conf_LatinName_ProtID_Dic_WithHuman,seq_file,file_name)#use protID_latin_name_dict to process the seq_file
				
				else: #failed cases are already handled in other functions
					pass


def orgranize_files(): #create folders to manage differernt types of files 
	current_dir = os.getcwd() + '/'

	for data_files in os.listdir(current_dir):
		if '.xml' in data_files and 'all_orthologues' in data_files:
			new_dir1 = current_dir + 'xml_files/'
			if not os.path.exists(new_dir1):
				os.makedirs(new_dir1)

		if 'latin_names' in data_files and '.fasta' in data_files:
			new_dir2 = current_dir + 'latinname_seqfiles/'
			if not os.path.exists(new_dir2):
				os.makedirs(new_dir2)

		if 'protIDs' in data_files and '.fasta' in data_files:
			new_dir3 = current_dir + 'proteinid_seqfiles/'
			if not os.path.exists(new_dir3):
				os.makedirs(new_dir3)

		if 'duplicated_speices' in data_files and '.txt' in data_files:
			new_dir4 = current_dir + 'dupspecies_record/'
			if not os.path.exists(new_dir4):
				os.makedirs(new_dir4)

		if '.txt' in data_files and 'failed' in data_files:
			new_dir5 = current_dir + 'fail_record/'
			if not os.path.exists(new_dir5):
				os.makedirs(new_dir5)

		if '.txt' in data_files and 'Failed_nohighor' in data_files:
			new_dir5 = current_dir + 'fail_record/'
			if not os.path.exists(new_dir5):
				os.makedirs(new_dir5)

		if '.txt' in data_files and 'missinginfo_' in data_files:
			new_dir5 = current_dir + 'fail_record/'
			if not os.path.exists(new_dir5):
				os.makedirs(new_dir5)

		if '.txt' in data_files and 'psudogene_' in data_files:
			new_dir5 = current_dir + 'fail_record/'
			if not os.path.exists(new_dir5):
				os.makedirs(new_dir5)

		if '.txt' in data_files and 'notenough_' in data_files:
			new_dir6 = current_dir + 'filtered_record/'
			if not os.path.exists(new_dir6):
				os.makedirs(new_dir6)

	for data_files in os.listdir(current_dir):
		if '.xml' in data_files and 'all_orthologues' in data_files:
			source = current_dir + data_files
			destination = new_dir1 + data_files
			shutil.move(source, destination)

		if 'latin_names' in data_files and '.fasta' in data_files:
			source = current_dir + data_files
			destination = new_dir2 + data_files
			shutil.move(source, destination)
					
		if 'protIDs' in data_files and '.fasta' in data_files:
			source = current_dir + data_files
			destination = new_dir3 + data_files
			shutil.move(source, destination)

		if 'duplicated_speices' in data_files and '.txt' in data_files:
			source = current_dir + data_files
			destination = new_dir4 + data_files
			shutil.move(source, destination)

		if '.txt' in data_files and 'failed' in data_files:
			source = current_dir + data_files
			destination = new_dir5 + data_files
			shutil.move(source, destination)

		if '.txt' in data_files and 'Failed_nohighor' in data_files:
			source = current_dir + data_files
			destination = new_dir5 + data_files
			shutil.move(source, destination)

		if '.txt' in data_files and 'missinginfo_' in data_files:
			source = current_dir + data_files
			destination = new_dir5 + data_files
			shutil.move(source, destination)

		if '.txt' in data_files and 'psudogene_' in data_files:
			source = current_dir + data_files
			destination = new_dir5 + data_files
			shutil.move(source, destination)

		if '.txt' in data_files and 'notenough_' in data_files:
			source = current_dir + data_files
			destination = new_dir6 + data_files
			shutil.move(source, destination)

#==================================================================================


if __name__ == '__main__':
	start_time = time.time()

	gene_file_name = sys.argv[2]
	path = sys.argv[1]
	check_or_not = sys.argv[3]

	stdout_name = gene_file_name.split('.txt')[0] + '_report.txt'
	symbol_file = open(gene_file_name,'r')
	symbol_lst = []
	for line in symbol_file:
		line = line.split('\n')[0] 
		symbol_lst.append(line) 
	print(str(len(symbol_lst)) + ' gene symbols queries')

	if check_or_not == 'y': #checked_dup already
		print('retrieving data...')
		print('output messages are redirected to ' + stdout_name + '...')
		sys.stdout = open(stdout_name, 'wr')
		for symbol_ in symbol_lst:
			main_operations(symbol_, path)

		print("FINISH ALL SYMBOL QUERY :)")
		print("Time Used: --- %s seconds ---" % (time.time() - start_time) + '\n')
		sys.stdout.close()
		orgranize_files()

	if check_or_not == 'n':#have not checked_duplication in gene file yet
	#-----------------------------------------------------------------
		#checking duplicated gene symbols in the txt file 
		checked = [] 
		duplicated = []
		for genes in symbol_lst:
			if genes not in checked:
				checked.append(genes)
			else:
				duplicated.append(genes)

		if(len(duplicated)>0):
			duplicated_file = open('duplicatedgenesymbols_'+ gene_file_name, 'w')
			duplicated_file.write(gene_file_name + '\n' + 'duplicated the following genes in gene_file\n'+str(duplicated))
			duplicated_file.close()
	#-----------------------------------------------------------------
		print('retrieving data...\n')
		print('output messages are redirect to ' + stdout_name + '...')
		sys.stdout = open(stdout_name, 'wr')
		for symbol_ in symbol_lst:
			main_operations(symbol_, path)
		print("FINISH ALL SYMBOL QUERY :)")
		print("Time Used: --- %s seconds ---" % (time.time() - start_time) + '\n')
		sys.stdout.close()
		orgranize_files()
