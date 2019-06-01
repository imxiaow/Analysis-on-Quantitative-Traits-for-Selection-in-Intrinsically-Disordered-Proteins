#!/usr/bin/python

#==================================================
#  Created on: Aug 2017
#      Author: xiao
#==================================================

#import libraries
import requests, sys
import json
import xml.etree.ElementTree as ET
import check_geneRefID_by_UniprotID as CKGENEID
import fetch_fetch_fetch as FFFROOT

#==================================================

#Defined variables 

server = "https://rest.ensembl.org" #REST Server
find_orthologues = "/homology/symbol/human/{}?" #REST API
symbol_look_up="/xrefs/name/human/{}" #REST API
look_up_iDs = "/xrefs/id/{}?" #REST API
get_sequence = "/sequence/id/{}?type=protein" #REST API

#Defined functions

def extract_high_conf_paired_geneRefID(root): #function use to extract high confidence ortholog in ortho file
	result = []
	high_conf_dict = {}
	gene_dict = {}

	for c in root.find('{http://orthoXML.org/2011/}groups'): #parsing xml orthologous file
		for c1 in c.findall('{http://orthoXML.org/2011/}property'):
			if c1.attrib['name'] == 'is_high_confidence':
				high_conf_dict[c.attrib['id']] = int(c1.attrib['value'])

		gene_dict[c.attrib['id']] = []
		for c1 in c.findall('{http://orthoXML.org/2011/}geneRef'):
			gene_dict[c.attrib['id']].append(c1.attrib['id'])

	for id_ in high_conf_dict:
		if high_conf_dict[id_] == 1:
			result.append(gene_dict[id_])

	if len(result) < 1:
		print("No get_high_conf_geneRefID avaliable in this gene file. Failed\n")

	return result

def get_high_conf_geneRefID(paired_geneRefID): #function use to get each species generef id
	GeneRef_Lst = []

	for pairs in paired_geneRefID:
		for numbers in pairs:
			if numbers not in GeneRef_Lst:
				GeneRef_Lst.append(numbers)

	return GeneRef_Lst

def re_extract_high_conf_geneRefID_H(human_geneID, root): 
	high_conf_paired_lst = extract_high_conf_paired_geneRefID(root)
	result = []

	for pairs in high_conf_paired_lst:
		if pairs[0] == human_geneID:
			result.append(pairs[1])

	result.append(human_geneID)

	print("Number of orthologues that is high_conf after re_extraction!!!!!" + str(len(result)))
	return result

def re_extract_Names_ProtID(root, geneRefId_Lst_H):
	names_ProtId_dict = {}
	geneRef_Id_homo = {}
	prot_Id_homo = []

	for c in root.findall('{http://orthoXML.org/2011/}species'):
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
							
	return names_ProtId_dict