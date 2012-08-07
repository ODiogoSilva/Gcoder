#!/usr/bin/python
# GCoder v1.0
# Author: Diogo N. Silva
# Last updated: 26/09/2011
# GCoder is a simple tool to convert indel events in an alignment of molecular sequences into a binary state block that can be interpreted by MrBayes (or any other software that allows for binary/restriction data). Currently, only Nexus format is supported as an output format.

#  Copyright 2012 Diogo N Silva <diogo@arch>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.


import argparse
import ElParsito
import re

parser = argparse.ArgumentParser(description="Concatenates DNA data matrices")
parser.add_argument("-if",dest="InputFormat",default="nexus",choices=["fasta","nexus","phylip"],help="Format of the input file (default is '%(default)s')")
parser.add_argument("-in",dest="infile",nargs=1,required=True,help="Input file. Only one allowed!")
parser.add_argument("-o",dest="outfile",default="Outfile",required=True,help="Name of the output file (default is '%(default)s')")

arg = parser.parse_args()

def dataset_creator (input_file,input_format,output_format):
	initial_storage,taxa_order = ElParsito.Taxa_gather(input_format,input_file)
	storage,part_list,sizes = ElParsito.Elparsito(input_format,initial_storage,input_file,output_format)
	return storage,taxa_order,part_list,sizes

def gap_listing (sequence,gap_symbol):
	""" Function that parses a sequence string and returns the position of indel events. The returned list is composed of tuples with the span of each indel """
	gap = "%s+" % (gap_symbol)
	span_regex = ""
	gap_list,seq_start = [],0
	while span_regex != None:
		span_regex = re.search(gap,sequence)
		if span_regex != None and seq_start == 0:
			gap_list.append(span_regex.span())
			sequence = sequence[span_regex.span()[1]+1:]
			seq_start = span_regex.span()[1]+1
		elif span_regex != None and seq_start != 0:
			gap_list.append((span_regex.span()[0]+seq_start,span_regex.span()[1]+seq_start))
			sequence = sequence[span_regex.span()[1]+1:]
			seq_start += span_regex.span()[1]+1
	return gap_list

def gap_binary_generator (sequence,gap_list):
	""" This function contains the algorithm to construct the binary state block for the indel events """
	inside_gap = "[a-z]?\-+[a-z]?"
	for cur_gap in gap_list:
		cur_gap_start,cur_gap_end = cur_gap
		if sequence[cur_gap_start:cur_gap_end] == "-"*(cur_gap_end - cur_gap_start) and sequence[cur_gap_start-1] != "-" and sequence[cur_gap_end] != "-":
			sequence += "1"
		elif sequence[cur_gap_start:cur_gap_end] == "-"*(cur_gap_end - cur_gap_start):
			if sequence[cur_gap_start-1] == "-" or sequence[cur_gap_end] == "-":
				sequence += "-"
		elif sequence[cur_gap_start:cur_gap_end] != "-"*(cur_gap_end - cur_gap_start):
			sequence += "0"
	return sequence

def multiSeq_gap_listing (storage,gap_symbol):
	complete_gap_list = []
	for taxa in storage:
		temp_list = gap_listing(storage[taxa],gap_symbol)
		unique_gap = [gap for gap in temp_list if gap not in complete_gap_list]
		complete_gap_list += unique_gap
	return complete_gap_list
	
def outfile_dump (storage,output_file,sizes,gap_list):
	out_file = open(arg.outfile+".nex","w")
	out_file.write("#NEXUS\n\nBegin data;\n\tdimensions ntax=%s nchar=%s;\n\tformat datatype=mixed(dna:1-%s,restriction:%s-%s) interleave=no gap=- missing=N;\n\tmatrix\n" % (len(storage),
						   sizes-1+len(gap_list),
						   sizes-1,
						   sizes,
						   len(gap_list)+sizes-1))
	for taxa in storage:
		out_file.write(taxa[:8].ljust(10)+" "+storage[taxa]+"\n")
	out_file.write("\t;\nend;\n\n[")
	for gap in gap_list:
		out_file.write("Indel event #%s: %s - %s\n" % (gap_list.index(gap)+1,gap[0],gap[1]))
	out_file.write("]")
	out_file.close()

def main():
	storage,taxa_order,part_list,sizes = dataset_creator(arg.infile,arg.InputFormat,"nexus")
	gap_list = multiSeq_gap_listing(storage,"-")
	for taxa in storage:
		storage[taxa] = gap_binary_generator (storage[taxa],gap_list)
	outfile_dump(storage,arg.outfile,sizes,gap_list)	

main()
