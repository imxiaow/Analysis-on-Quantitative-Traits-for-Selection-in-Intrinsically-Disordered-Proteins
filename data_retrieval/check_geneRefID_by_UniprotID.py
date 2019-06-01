#!/usr/bin/python

#==================================================
#  Created on: Aug 2017
#      Author: xiao
#==================================================

#import libraries
import requests, sys
import json
import time

#==================================================

#Defined variables 

server = "https://rest.ensembl.org" #REST Server
look_up_iDs = "/xrefs/id/{}?" #REST API
symbol_look_up="/xrefs/name/human/{}" #REST API

#Defined functions

def get_info_UNPIT(server,request): #function use to get Uniprot ID from Ensembl GeneRefID
	r = requests.get(server+request, headers={ "Content-Type" : "application/json"})

	if not r.ok:
		error_code = r.status_code
  		if error_code == 429: #too many request, reached rate limit
  			print("wait..\n")
  			time.sleep(0.5)
  			return get_info_UNPIT(server,request)

  		else:
  			print("error in get_info_UNIPIT\n")
  			return -1

  	else:
		decoded =r.json()
		infos = {}

		for item in decoded:
			if 'db_display_name' in item and item['db_display_name'] == 'UniProtKB Gene Name':
				infos[str(item['primary_id'])]=request.strip('/xrefsid?')

		print(infos)
		return infos

def get_symbol_UNIPT_info(server,request): #function use to get Uniprot ID from Gene Symbol 
	r = requests.get(server+request, headers={ "Content-Type" : "application/json"})

	if not r.ok:
		error_code = r.status_code
  		if error_code == 429:
  			print("wait..\n")
  			time.sleep(0.5)
  			return get_symbol_UNIPT_info(server,request)

  		else:
  			print("error in get_symbol_UNIPT_info\n")
  			return -2

  	else:
		decoded =r.json()
		infos = {}

		for item in decoded:
			if 'db_display_name' in item and item['db_display_name'] == 'UniProtKB Gene Name': #parsing json orthologous file
				infos[str(item['primary_id'])]=request.strip('/xrefs/name/human/')

		print(infos)
		return infos

def get_info_HGNC(server, request): # function use to get HGNC symbol information for ids
	r = requests.get(server+request, headers={ "Content-Type" : "application/json"})

	if not r.ok:
		error_code = r.status_code
  		if error_code == 429:
  			print("wait..\n")
  			time.sleep(0.5)
  			return get_info_HGNC(server, request)

  		else:
  			print("error in get_info_HGNC\n")
 			return -3

	else:
		decoded =r.json()
		HGNC_Id = {}

		for item in decoded:
			if 'dbname' in item and item['dbname'] == 'HGNC': #parsing json orthologous file
				HGNC_Id[str(item['primary_id'])]=request.strip('/xrefsid?')

		print(HGNC_Id)
		return HGNC_Id

def get_symbol_HGNC(server,request): #function use to get HGNC symbol from gene name
	r = requests.get(server+request, headers={ "Content-Type" : "application/json"})

	if not r.ok:
		error_code = r.status_code
  		if error_code == 429:
  			print("wait..\n")
  			time.sleep(0.5)
  			return get_symbol_HGNC(server,request)

  		else:
  			print("error in get_symbol_HGNC\n")
 			return -4

	else:
		decoded =r.json()
		HGNC_symbol = {}

		for item in decoded:
			if 'dbname' in item and item['dbname'] == 'HGNC': #parsing json orthologous file
				HGNC_symbol[str(item['primary_id'])]=request.strip('/xrefs/name/human/')

		print(HGNC_symbol)
		return HGNC_symbol

def UNIPT_info_compare(symbol, Id_, server,symbol_look_up, look_up_iDs): #function comparing uniprot information
	check_get_symbol_HGNC = get_symbol_HGNC(server,symbol_look_up.format(symbol))

	if check_get_symbol_HGNC == -4:
		write_error_getsymbolHGNC_filename = 'error_getsymbolHGNC_' + symbol + '.txt'
		write_error_getsymbolHGNC_file = open(write_error_getsymbolHGNC_filename,'wr')
		write_error_getsymbolHGNC_file.write(symbol+"\nerror in get_symbol_HGNC\n")
		write_error_getsymbolHGNC_file.close()
		return -41

	else:
		symbol_HGNC_info = check_get_symbol_HGNC
		if symbol_HGNC_info == {}:
			return -1 # no HGNC info for query symbol

		elif len(symbol_HGNC_info) >1:
			return -5

		else: #exactly one HGNC for symbol
			Ids_HGNC_info = []
			for i in Id_:
				check_get_ID_HGNC = get_info_HGNC(server, look_up_iDs.format(i))
				if check_get_ID_HGNC == -3:
					write_error_getinfoHGNC_filename = 'error_getinfoHGNC_' + symbol + '.txt'
					write_error_getinfoHGNC_file = open(write_error_getinfoHGNC_filename,'wr')
					write_error_getinfoHGNC_file.write(symbol+'\n'+i+"\nerror in get_info_HGNC\n"+'the id that need distinguished')
					write_error_getinfoHGNC_file.close()
					return -42

				if check_get_ID_HGNC != {}:
					Ids_HGNC_info.append(check_get_ID_HGNC)
				
		if len(Ids_HGNC_info) < 1: # no HGNC info for all IDs
			return -2

		else:
			result_real_geneID_HGNC = []
			for symbol_HGNC in symbol_HGNC_info:
				for Ids_dict in Ids_HGNC_info:
					for Id_HGNC in Ids_dict:
						if Id_HGNC == symbol_HGNC:
							if Ids_dict[Id_HGNC] not in result_real_geneID_HGNC:
								result_real_geneID_HGNC.append(Ids_dict[Id_HGNC])

			#now check comparison result
			if len(result_real_geneID_HGNC) == 1:
				return result_real_geneID_HGNC

			elif len(result_real_geneID_HGNC) == 0:
				return -3 # no matched HGNC between Ids and query symbol 

			else:
				return -4 # can not distinguish bwtween Ids and query symbol by HGNC
				

def get_real_symbol_geneID(symbol, Id_, server,symbol_look_up, look_up_iDs): #function handle all the cases would appear in gene names

	check_get_symbol_UNIPT_info = get_symbol_UNIPT_info(server,symbol_look_up.format(symbol))

	if check_get_symbol_UNIPT_info == -2:
		write_error_getsymbolUNIPT_filename = 'error_getsymbolUNIPT_' + symbol + '.txt'
		write_error_getsymbolUNIPT_file = open(write_error_getsymbolUNIPT_filename,'wr')
		write_error_getsymbolUNIPT_file.write(symbol+"\nerror in get_symbol_UNIPT_info\n")
		write_error_getsymbolUNIPT_file.close()
		return -2

	else:
		symbol_UNIPT_info = check_get_symbol_UNIPT_info
		if symbol_UNIPT_info == {}: #no Uniprot info for symbol
			check_UNIPT_info_compare = UNIPT_info_compare(symbol, Id_, server,symbol_look_up, look_up_iDs)
			if check_UNIPT_info_compare == -41:
				return 1234

			if check_UNIPT_info_compare == -42:
				return 1234

			if check_UNIPT_info_compare == -1:
				write_failfilename= 'failed_' + symbol + '.txt'
				write_failfile = open(write_failfilename,'wr')
				write_failfile.write(symbol + '\n' + "No Uniprot for query symbol\nNo HGNC for query symbol\n")
				write_failfile.close()
				return 1234

			elif check_UNIPT_info_compare == -2:
				write_failfilename= 'failed_' + symbol + '.txt'
				write_failfile = open(write_failfilename,'wr')
				write_failfile.write(symbol + '\n' + "No Uniprot for query symbol\nNo IDshave HGNC\n")
				write_failfile.close()
				return 1234

			elif check_UNIPT_info_compare == -3:
				write_failfilename= 'failed_' + symbol + '.txt'
				write_failfile = open(write_failfilename,'wr')
				write_failfile.write(symbol + '\n' + "No Uniprot for query symbol\nNo matched HGNC\n")
				write_failfile.close()
				return 1234

			elif check_UNIPT_info_compare == -4:
				write_failfilename= 'failed_' + symbol + '.txt'
				write_failfile = open(write_failfilename,'wr')
				write_failfile.write(symbol + '\n' + "No Uniprot for query symbol\nMore than one matched HGNC\ncan't distinguish\n")
				write_failfile.close()
				return 1234

			elif check_UNIPT_info_compare == -5:
				write_failfilename= 'failed_' + symbol + '.txt'
				write_failfile = open(write_failfilename,'wr')
				write_failfile.write(symbol + '\n' + "No Uniprot for query symbol\nMore than one HGNC for query symbol\ncant compare")
				write_failfile.close()
				return 1234

			else:
				real_symbol_geneID = check_UNIPT_info_compare
				return real_symbol_geneID

		else: # symbol have Uniprot ID
			if len(symbol_UNIPT_info) > 1: 
				check_UNIPT_info_compare = UNIPT_info_compare(symbol, Id_, server,symbol_look_up, look_up_iDs)

				if check_UNIPT_info_compare == -41:
					return 5678

				if check_UNIPT_info_compare == -42:
					return 5678

				if check_UNIPT_info_compare == -1:
					write_failfilename= 'failed_' + symbol + '.txt'
					write_failfile = open(write_failfilename,'wr')
					write_failfile.write(symbol + '\n' + "More than one UNPROT for this symbol\nNo HGNC for symbol\n")
					write_failfile.close()
					return 5678

				elif check_UNIPT_info_compare == -2:
					write_failfilename= 'failed_' + symbol + '.txt'
					write_failfile = open(write_failfilename,'wr')
					write_failfile.write(symbol + '\n' + "More than one UNPROT for this symbol\nNo IDshave HGNC\n")
					write_failfile.close()
					return 5678

				elif check_UNIPT_info_compare == -3:
					write_failfilename= 'failed_' + symbol + '.txt'
					write_failfile = open(write_failfilename,'wr')
					write_failfile.write(symbol + '\n' + "More than one UNPROT for this symbol\nNo matched HGNC\n")
					write_failfile.close()
					return 5678

				elif check_UNIPT_info_compare == -4:
					write_failfilename= 'failed_' + symbol + '.txt'
					write_failfile = open(write_failfilename,'wr')
					write_failfile.write(symbol + '\n' + "More than one UNPROT for this symbol\nMore than one matched HGNC\ncan't distinguish\n")
					write_failfile.close()
					return 5678

				elif check_UNIPT_info_compare == -5:
					write_failfilename= 'failed_' + symbol + '.txt'
					write_failfile = open(write_failfilename,'wr')
					write_failfile.write(symbol + '\n' + "More than one UNPROT for this symbol\nMore than one HGNC for query symbol\ncant compare")
					write_failfile.close()
					return 5678

				else:
					real_symbol_geneID = check_UNIPT_info_compare
					return real_symbol_geneID

			else: # symbol have exactly one Uniprot ID
				Ids_UNIPT_info = []
				for i in Id_:
					check_get_ID_UNIPT = get_info_UNPIT(server, look_up_iDs.format(i))
					if check_get_ID_UNIPT == -3:
						write_error_getinfoUNIPT_filename = 'error_getinfoUNPIT' + symbol + '.txt'
						write_error_getinfoUNIPT_file = open(write_error_getinfoUNIPT_filename,'wr')
						write_error_getinfoUNIPT_file.write(symbol+'\n'+i+"\nerror in get_info_UNPIT\n"+'the id that need distinguished')
						write_error_getinfoUNIPT_file.close()
						return 9000

					if check_get_ID_UNIPT != {}:
						Ids_UNIPT_info.append(check_get_ID_UNIPT)
	
				if len(Ids_UNIPT_info) < 1:# no UNIPROT info for all IDs
					check_UNIPT_info_compare = UNIPT_info_compare(symbol, Id_, server,symbol_look_up, look_up_iDs)
					if check_UNIPT_info_compare == -41:
						return 9000

					if check_UNIPT_info_compare == -42:
						return 9000

					if check_UNIPT_info_compare == -1:
						write_failfilename= 'failed_' + symbol + '.txt'
						write_failfile = open(write_failfilename,'wr')
						write_failfile.write(symbol + '\n' + "There are exactly one UNIPROT for this symbol\nno UNIPROT info for all IDs\nNo HGNC for symbol\n")
						write_failfile.close()
						return 9000

					elif check_UNIPT_info_compare == -2:
						write_failfilename= 'failed_' + symbol + '.txt'
						write_failfile = open(write_failfilename,'wr')
						write_failfile.write(symbol + '\n' + "There are exactly one UNIPROT for this symbol\nno UNIPROT info for all IDs\nNo IDshave HGNC\n")
						write_failfile.close()
						return 9000

					elif check_UNIPT_info_compare == -3:
						write_failfilename= 'failed_' + symbol + '.txt'
						write_failfile = open(write_failfilename,'wr')
						write_failfile.write(symbol + '\n' + "There are exactly one UNIPROT for this symbol\nno UNIPROT info for all IDs\nNo matched HGNC\n")
						write_failfile.close()
						return 9000

					elif check_UNIPT_info_compare == -4:
						write_failfilename= 'failed_' + symbol + '.txt'
						write_failfile = open(write_failfilename,'wr')
						write_failfile.write(symbol + '\n' + "There are exactly one UNIPROT for this symbol\nno UNIPROT info for all IDs\nMore than one matched HGNC\ncan't distinguish")
						write_failfile.close()
						return 9000

					elif check_UNIPT_info_compare == -5:
						write_failfilename= 'failed_' + symbol + '.txt'
						write_failfile = open(write_failfilename,'wr')
						write_failfile.write(symbol + '\n' + "There are exactly one UNIPROT for this symbol\nno UNIPROT info for all IDs\nMore than one HGNC for query symbol\ncant compare")
						write_failfile.close()
						return 9000

					else:
						real_symbol_geneID = check_UNIPT_info_compare
						return real_symbol_geneID

				else: # Ids have Uniprot ID
					result_real_geneID_UNIPT = []
					for symbol_UNIPT in symbol_UNIPT_info:
						for Ids_dict in Ids_UNIPT_info:
							for Id_UNIPT in Ids_dict:
								if Id_UNIPT == symbol_UNIPT:
									if Ids_dict[Id_UNIPT] not in result_real_geneID_UNIPT:
										result_real_geneID_UNIPT.append(Ids_dict[Id_UNIPT])

					if len(result_real_geneID_UNIPT) == 1:
						return result_real_geneID_UNIPT

					elif len(result_real_geneID_UNIPT) == 0: #no matach
						check_UNIPT_info_compare = UNIPT_info_compare(symbol, Id_, server,symbol_look_up, look_up_iDs)
						if check_UNIPT_info_compare == -41:
							return 91234

						if check_UNIPT_info_compare == -42:
							return 91234

						if check_UNIPT_info_compare == -1:
							write_failfilename= 'failed_' + symbol + '.txt'
							write_failfile = open(write_failfilename,'wr')
							write_failfile.write(symbol + '\n' + "There are exactly one UNIPROT for this symbol\nNo matched UNPROT\nNo HGNC for symbol\n")
							write_failfile.close()
							return 91234

						elif check_UNIPT_info_compare == -2:
							write_failfilename= 'failed_' + symbol + '.txt'
							write_failfile = open(write_failfilename,'wr')
							write_failfile.write(symbol + '\n' + "There are exactly one UNIPROT for this symbol\nNo matched UNPROT\nNo IDshave HGNC\n")
							write_failfile.close()
							return 91234

						elif check_UNIPT_info_compare == -3:
							write_failfilename= 'failed_' + symbol + '.txt'
							write_failfile = open(write_failfilename,'wr')
							write_failfile.write(symbol + '\n' + "There are exactly one UNIPROT for this symbol\nNo matched UNPROT\nNo matched HGNC\n")
							write_failfile.close()
							return 91234

						elif check_UNIPT_info_compare == -4:
							write_failfilename= 'failed_' + symbol + '.txt'
							write_failfile = open(write_failfilename,'wr')
							write_failfile.write(symbol + '\n' + "There are exactly one UNIPROT for this symbol\nNo matched UNPROT\nMore than one matched HGNC\ncan't distinguish")
							write_failfile.close()
							return 91234

						elif check_UNIPT_info_compare == -5:
							write_failfilename= 'failed_' + symbol + '.txt'
							write_failfile = open(write_failfilename,'wr')
							write_failfile.write(symbol + '\n' + "There are exactly one UNIPROT for this symbol\nNo matched UNPROT\nMore than one HGNC for query symbol\ncant compare")
							write_failfile.close()
							return 91234

						else:
							real_symbol_geneID = check_UNIPT_info_compare
							return real_symbol_geneID

					else: # more than one match
						check_UNIPT_info_compare = UNIPT_info_compare(symbol, Id_, server,symbol_look_up, look_up_iDs)
						if check_UNIPT_info_compare == -41:
							write_failfilename= 'failed_' + symbol + '.txt'
							write_failfile = open(write_failfilename,'wr')
							write_failfile.write(symbol + '\n' + "There are exactly one UNIPROT for this symbol\nMore than one matched UNPROT\n")
							write_failfile.close()
							return 95678

						if check_UNIPT_info_compare == -42:
							write_failfilename= 'failed_' + symbol + '.txt'
							write_failfile = open(write_failfilename,'wr')
							write_failfile.write(symbol + '\n' + "There are exactly one UNIPROT for this symbol\nMore than one matched UNPROT\n")
							write_failfile.close()
							return 95678

						if check_UNIPT_info_compare == -1:
							write_failfilename= 'failed_' + symbol + '.txt'
							write_failfile = open(write_failfilename,'wr')
							write_failfile.write(symbol + '\n' + "There are exactly one UNIPROT for this symbol\nMore than one matched UNPROT\nNo HGNC for symbol\n")
							write_failfile.close()
							return 95678

						elif check_UNIPT_info_compare == -2:
							write_failfilename= 'failed_' + symbol + '.txt'
							write_failfile = open(write_failfilename,'wr')
							write_failfile.write(symbol + '\n' + "There are exactly one UNIPROT for this symbol\nMore than one matched UNPROT\nNo IDshave HGNC\n")
							write_failfile.close()
							return 95678

						elif check_UNIPT_info_compare == -3:
							write_failfilename= 'failed_' + symbol + '.txt'
							write_failfile = open(write_failfilename,'wr')
							write_failfile.write(symbol + '\n' + "There are exactly one UNIPROT for this symbol\nMore than one matched UNPROT\nNo matched HGNC\n")
							write_failfile.close()
							return 95678

						elif check_UNIPT_info_compare == -4:
							write_failfilename= 'failed_' + symbol + '.txt'
							write_failfile = open(write_failfilename,'wr')
							write_failfile.write(symbol + '\n' + "There are exactly one UNIPROT for this symbol\nMore than one matched UNPROT\nMore than one matched HGNC\ncan't distinguish")
							write_failfile.close()
							return 95678

						elif check_UNIPT_info_compare == -5:
							write_failfilename= 'failed_' + symbol + '.txt'
							write_failfile = open(write_failfilename,'wr')
							write_failfile.write(symbol + '\n' + "There are exactly one UNIPROT for this symbol\nMore than one matched UNPROT\nMore than one HGNC for query symbol\ncant compare")
							write_failfile.close()
							return 95678

						else:
							real_symbol_geneID = check_UNIPT_info_compare
							return real_symbol_geneID