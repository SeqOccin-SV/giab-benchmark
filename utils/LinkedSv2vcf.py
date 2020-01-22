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
0	 1				2	 3			4		  5					6			7		  8		9	 10
Chr1 Start1			Chr2 Start2		#SplitMol #DiscordantReads	orientation	Haplotype score pass filter
'''

''' VCF FORMAT
For LR
chr pos ID		? type		QUAL FILTER		INFOS	GT 0|1

For GIAB

For pipe ONF
Chrom pos ID	REF alt 	Qual filter		INFOS	Format Sample
'''

# Basically, everyone is doing smth different in bedpe, vcf is close to being standardize

# exemple header for VCF
headerString = '''##fileformat=VCFv4.2
##source=LinkedSV_Conversion
##FILTER=<ID=PASS,Description="Passed filter from software">
##FILTER=<ID=FAIL,Description="Failed filter form software">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of SV:DEL=Deletion, CON=Contraction, INS=Insertion, DUP=Duplication, INV=Inversion">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GTcons1,Number=1,Type=String,Description="Consensus Genotype using the GT from svviz2 rather than ref and alt allele counts, which is sometimes inaccurate for large variants">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE'''

print(headerString)
# Lines list. Needed to output the lines in numerical order at the end
lines = []

# Reading LinkedSV bedpe file
with open(inputname, 'rt') as fh:
	for line in fh:
		# Skip comment
		if (line[0] == '#'):
			continue
		# Attempt to skip empty line
		#if (line[0] == ''):
		#	continue
		line_content = line.rstrip().split('\t')
		
		# Skipping empty line
		if (len(line_content) < 9):
			continue
			
		#info_content = line_content[11].split(';')
		#infos = dict()
		#for info in info_content:
		#	key, data = info.split('=')
		#	infos[key] = data
		#print(json.dumps(infos, indent=1))
		#break
		
		#svlen = str(abs(int(line_content[4])-int(line_content[1])))
		
		output_infos = ';'.join( ('END='+line_content[4], 'SVTYPE='+line_content[6], 'SVLEN='+line_content[8] ) )
		output_line = '\t'.join((line_content[0],line_content[1],line_content[7],'.','<'+line_content[6]+'>',line_content[9],line_content[10],output_infos,'.','.'))
		
		lines.append(output_line)
		# break

# Final alphanumeric sort and print
output = sort_human(lines)
for line in output:
	print(line)
