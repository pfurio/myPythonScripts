#!/usr/bin/python
# -*- coding: utf-8 -*- 

###########   2015-07-23 pfurio@cipf.es             ###########
# Input: miRNA.dat input in EMBL format (http://www.mirbase.org/ftp.shtml) and the species to be extracted
# Returns: Mature - Precursor miRNA - Position within the precursor
# Execution example: ./extract_mirnas.py miRNA.dat hsa hsa_miRNAs_tab.txt

import sys
import re

mirna   = open(sys.argv[1],'r')
species = sys.argv[2]
miout   = open(sys.argv[3],'w')

flag    = False
matures = ""
positions = ""
for line in mirna:
	if line[:2] == "ID" and re.search(species,line) is not None:
		flag = True
		precursor = line.split()[1]
		matures = ""
		positions = ""
	elif flag is True and line[:2] == "FT" and re.search("miRNA",line) is not None:
		positions = line.split()[2]
	elif flag is True and line[:2] == "FT" and re.search("product",line) is not None:
		matures = line.split("\"")[1]
		miout.write(precursor + "\t" + matures + "\t" + '\t'.join(positions.split("..")) + "\n")	
	elif flag is True and line[:2] == "SQ":
		flag = False
		
mirna.close()
miout.close()

