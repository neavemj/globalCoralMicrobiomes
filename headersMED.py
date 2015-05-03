#!/usr/bin/env python

## Script to take mothur fasta and groups file 
## and incorporate information suitable for minimum entropy decomposition

# Matthew J. Neave 23.12.14 **matthewjneave1@gmail.com**

import sys
import argparse


# parse arguments from the command line

parser = argparse.ArgumentParser("Convert mothur fasta and groups files to MED compatible format")

parser.add_argument('mothur_fasta_file', type = argparse.FileType("r"),
		nargs = "?", help = "fasta file from mothur")
parser.add_argument('mothur_groups_file', type = argparse.FileType("r"),
		nargs = "?", help = "groups file from mothur")
parser.add_argument('output_file_name', type = argparse.FileType("w"), 
		nargs = "?", help = "output file name")

args = parser.parse_args()

out_handle = args.output_file_name

groupsDict = {}

for line in args.mothur_groups_file:
	cols = line.split()
	seq = cols[0]
	group = cols[1]
	groupsDict[seq] = group


missingHeader = 0

for line in args.mothur_fasta_file:
	line = line.strip()
	if line.startswith('>'):
		header = line.lstrip('>')
		if header in groupsDict:
			newHeader = header.replace('_', ',')
			out_handle.write('>' + groupsDict[header] + '_' + newHeader + '\n')
		else:
			missingHeader +=1
	else:
		out_handle.write(line + '\n')


if missingHeader > 0:
	print 'headers not found: ', missingHeader
else:
	print 'all headers successfully converted'


