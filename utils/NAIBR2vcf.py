#!/usr/bin/env python
# -*-coding:utf-8 -*

import os.path
import sys
import json
import re

# a function to sort in alphanumeric
def sort_human(l):
    convert = lambda text: float(text) if text.isdigit() else text
    alphanum = lambda key: [convert(c) for c in re.split('([-+]?[0-9]*\.?[0-9]*)', key)]
    l.sort(key=alphanum)
    return l


try:
	inputname = sys.argv[1]
	os.path.isfile(inputname)
except IndexError:
	print('You need to give an input bedpe file as argument')
except IOError:
	print('Unable to open user input file')

''' BEDPE FORMAT
For LR,
0	1	2			3	4	5			6	7		8		9		10		11
chr1 start1 end1	chr2 start2 end2	ID Score	Strand1 Strand2 filters	INFO

For LinkedSV
0	1	2			3	4	5			6	7	8		9		10	11
chr1 start1 end1	chr2 start2 end2	Type ID Length	Score filter INFO=Barcode

For NAIBR
0	 1				2	 3			4		  5					6			7		  8		9
Chr1 Start1			Chr2 Start2		#SplitMol #DiscordantReads	orientation	Haplotype score PASS_filter
'''

''' VCF FORMAT
For LR
chr pos ID		? type		QUAL FILTER		INFOS	GT 0|1

For GIAB

For pipe ONF
Chrom pos ID	REF alt 	Qual filter		INFOS	Format Sample
'''

# Basically, everyone is doing smth different in bedpe, vcf is close to being standardize

header_string = '''##fileformat=VCFv4.2
##source=NAIBR_bedpe_conversion
##FILTER=<ID=PASS,Description="Passed the software filter">
##FILTER=<ID=FAIL,Description="Failed the software filter">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of SV:DEL=Deletion, CON=Contraction, INS=Insertion, DUP=Duplication, INV=Inversion">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GTcons1,Number=1,Type=String,Description="Consensus Genotype using the GT from svviz2 rather than ref and alt allele counts, which is sometimes inaccurate for large variants">
##FORMAT=<ID=PB_GT,Number=1,Type=String,Description="Genotype predicted by svviz from PacBio">
##FORMAT=<ID=PB_REF,Number=1,Type=Integer,Description="Number of PacBio reads supporting the REF allele as predicted by svviz">
##FORMAT=<ID=PB_ALT,Number=1,Type=Integer,Description="Number of PacBio reads supporting the ALT allele as predicted by svviz">
##FORMAT=<ID=PBHP_GT,Number=1,Type=String,Description="Genotype predicted by svviz from PacBio haplotype-separated by 10X">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE'''

# Lines list. Needed to output the lines in numerical order at the end
lines = []

# Reading NAIBR bedpe file and parsing each line
with open(inputname, 'rt') as fh:
	i=0
	for line in fh:
		# Skip header
		if i == 0:
			i += 1
			print(header_string)
			continue
		line_content = line.rstrip().split('\t')
		
		# Skip infos, useless for NAIBR
		#print(line_content)
		#print(json.dumps(line_content, indent=1))
		
		#info_content = line_content[11].split(';')
		#infos = dict()
		#for info in info_content:
		#	key, data = info.split('=')
		#	infos[key] = data
		#print(json.dumps(infos, indent=1))
		#break
		
		# Retrieving sv type information
		orientation = line_content[6]
		svtype = ''
		if (orientation == '+-' or orientation == '-+' ):
			svtype = 'INV'
		elif (orientation == '--'):
			svtype = 'DEL'
		elif (orientation == '++'):
			svtype = 'INS'
			
		# Construction of some fields like ID and INFO
		svlen = str(abs(int(line_content[3])-int(line_content[1])))
		name = 'NAIBR_'+str(i)
		output_infos = ';'.join( ('END='+line_content[3], 'SVTYPE='+svtype, 'SVLEN='+svlen ) )
		# Final line construction and appending to list
		output_line = '\t'.join((line_content[0],line_content[1],name,'.','<'+svtype+'>',line_content[8],line_content[9],output_infos,'.','.'))
		lines.append(output_line)
		i += 1
		# break

# Final alphanumeric sort and print
output = sort_human(lines)
for line in output:
	print(line)
