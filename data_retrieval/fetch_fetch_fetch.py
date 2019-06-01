#!/usr/bin/python

#==================================================
#  Created on: Aug 2017
#      Author: xiao
#==================================================

#import libraries
import requests, sys
import json
import os
import time
import xml.etree.ElementTree as ET

#==================================================

#Defined variables 

server = "https://rest.ensembl.org" #REST Server
find_orthologues = "/homology/symbol/human/{}?" #REST API
symbol_look_up="/xrefs/name/human/{}" #REST API
look_up_iDs = "/xrefs/id/{}?" #REST API
get_sequence = "/sequence/id/{}?type=protein" #REST API
 
#Defined functions 

def fetch_orthologues(server, request, write_file): #function use REST API to get ortholog information
	r = requests.get(server+request, headers={ "Content-Type" : "text/x-orthoxml+xml"})
	if not r.ok:
  		error_code = r.status_code
  		if error_code == 400:
  			print("error 400 fetch_ortholog\n")
  			return -1

  		if error_code == 429: #too many request, reached rate limit
  			time.sleep(0.5)
  			print("waiting\n") 
  			return fetch_orthologues(server, request, write_file)

  		else:
  			print("error in fetch_orthologues request\n")
  			return -2

  	else:
  		write_file.write(r.text)

def write_process_file_getroot(command, server): #function use to process xml orthologue file, check data before proceed
	write_file_name = command + '_all_orthologues.xml'
	write_all_ortho_file = open(write_file_name,'wr')
	check = fetch_orthologues(server, find_orthologues.format(command), write_all_ortho_file)
	if check == -1:
		print("No Information found by this symbol: " + command+ '\ncheck '+ command + ' on Ensembl web.')
		write_miss_infofilename = 'missinginfo_' + command + '.txt'
		write_missfile = open(write_miss_infofilename,'wr')
		write_missfile.write(command + '\n' + "No Information found by this symbol, check this symbol on Ensembl web." +'\n')
		write_missfile.close()
		return None

	if check == -2:
		write_error_filename = 'error_fetchortho_' + command + '.txt'
		write_erorfile = open(write_error_filename, 'wr')
		write_erorfile.write(command + '\n' + "error in fetch_orthologues request\n")
		write_erorfile.close()
		return None

	if check == None:
		write_all_ortho_file.close()
		tree = ET.parse(write_file_name)
		root = tree.getroot()
		return root

def fetch_sequence(server, request): #function use REST API to get sequence 
	r = requests.get(server+request, headers={ "Content-Type" : "text/x-fasta"})
	if not r.ok:
  		error_code = r.status_code
  		if error_code == 429: #too many request, reached rate limit
  			time.sleep(0.5)
  			return fetch_sequence(server, request)

  		else:
  			print("error in fetch_sequence request\n")
  			return -1

  	else:
		return r.text

def process_name_protID_dict_to_lst(name_protID_dict): #helper function to get all protein ID in dictionary
	lst=[]
	for name in name_protID_dict:
		for item in name_protID_dict[name]:
			lst.append(item)

	return lst

def write_sequence_file(name_protID_dict,command): # function to write sequence file
	write_file_name = 'all_seqs_with_protIDs$'+command + '.fasta'
	write_file = open(write_file_name,'wr')
	protID_lst = process_name_protID_dict_to_lst(name_protID_dict)
	#print(str(len(protID_lst))+"\n")
	check = 0
	for protID in protID_lst:
		check_seqinfo = fetch_sequence(server,get_sequence.format(protID))
		if check_seqinfo == -1:
			write_missseq_filename = 'error_fetchseq' + command + '_' + str(protID) + '.txt'
			write_missseqfile = open(write_missseq_filename,'wr')
			write_missseqfile.write(command + '\n' + "error in fetch_sequence request\n" + protID)
			write_missseqfile.close()

		elif check_seqinfo != -1:
			sequence = check_seqinfo
			write_file.write(sequence)
			check +=1

	print("checking number of sequences writing to file: "+ str(check) +'\n')
	write_file.close()
	print("Finish writing sequence file: " + write_file_name)

def protID_latin_name_dict(name_protID_dict): # helper function to exchange key and values in dictionary
	protID_name_dict = {}
	for name in name_protID_dict:
		for protID in name_protID_dict[name]:
			protID_name_dict[protID] = name

	return protID_name_dict

def process_seq_file_to_latinName_wt(name_protID_dict, seq_file, seq_filename): #function use to convert protein id in sequence file into latin name sequence file
	file_name = seq_filename.strip('.fasta').split('$')[1]
	write_file_name = 'all_seqs_with_latin_names$' + file_name +'.fasta'
	write_file = open(write_file_name,'wr')
	dict_protID_latinname = protID_latin_name_dict(name_protID_dict)

	for line in seq_file:
		if line.startswith('>'):
			line = line.strip('>').split('\n')[0]
			if line in dict_protID_latinname:
				latin_name = dict_protID_latinname[line]
				new_line = '>' + latin_name + '\n'
				write_file.write(new_line)

			else:
				write_file.write("KeyError!!!!!!!!!!!!!!\n")

		else:
			write_file.write(line)
			
	write_file.close()
	print('Finish converting sequence file to latin name: '+ write_file_name)	