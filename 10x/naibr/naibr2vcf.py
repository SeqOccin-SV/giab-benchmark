#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
  Converting NAIBR bedpe to vcf format

"""

import os
import argparse
import datetime

from collections import defaultdict
import collections
from pysam import VariantFile
from pybedtools import BedTool


_Contig = collections.namedtuple('Contig', ['name', 'length'])


def get_contigs(reference):
    """
    Retrieve the contigs (chromosomes) from the genome fasta file
    :param reference: the genome fasta file
    :type reference: fasta filefile
    :return: list of contigs in the format specified by PyCvf
    :rtype: list
    """
    contigs = []
    reference_fai = str(reference) + ".fai"
    if reference is not None and os.path.isfile(reference_fai):
        with open(reference_fai) as fai_file:
            for line in fai_file:
                line_items = line.strip().split("\t")
                name, length = line_items[0:2]
                name = name.split(" ")[0]
                contigs.append(_Contig(name, int(length)))
    return contigs


def vcf_from_scratch(output_file, contigs):
    """
    A tiny vcf file whith a minimal header with contig information
    """
    date = datetime.datetime.now()
    tmp_vcf = output_file + "_tmp"
    with open(tmp_vcf, "w") as fout:
        fout.write("##fileformat=VCFv4.2\n")
        fout.write("##fileDate=%s\n" % date.strftime("%Y%M%d"))
        for ctg in contigs:
            fout.write("##contig=<ID=%s,length=%d>\n" % (ctg.name, ctg.length))
        fout.write("#%s\n" % "\t".join(["CHROM", "POS", "ID", "REF", "ALT",
                                        "QUAL", "FILTER", "INFO", "FORMAT"]))
    return tmp_vcf


def reformat_bedpe2vcfrecords(bedpefile, vcf_header):
    """
        Parsing bedpe records to construct the records
        Genotype is set to unknown
    """
    # naibr format is not recognize by BedTool lib
    #naibr_variants = BedTool(bedpefile)
    
    # Bedpe fields : https://github.com/raphael-group/NAIBR
    keys = ['chrom1', 'break1',
            'chrom2', 'break2',
            'NSplitMol', 'NDiscordant', 'orientation',
            'haplotype', 'score', 'filter']

    records = []
    infos = defaultdict()
    i = 0
    
    with open(bedpefile) as naibr_variants:
        for row in naibr_variants:
            # skipping first row
            if (i==0):
                i += 1
                continue
            variant = dict(zip(keys, row.rstrip().split('\t')))
            samples = [{'GT': (None, None)}]
            '''
            From NAIBR github: https://github.com/raphael-group/NAIBR/issues/11
            "The orientation field is actually key to distinguishing between SV types."
            "In the reference genome segments are aligned head to tail, so '+-' orientation is the reference orientation (which includes deletions). '--' and '++' indicate an novel adjacency created by a inversion and '-+' indicates a novel adjacency created by a tandem duplication."
            +- -> Deletion (and insertion?)
            -- | ++ -> Inversion
            -+ -> Tandem Duplication
            '''
            # Filtering to conserve only deletion
            print(str(variant['orientation']))
            if (str(variant['orientation']) == '+-'):
                print('success')
                infos.clear()
                infos['SVTYPE'] = 'DEL'
                infos['SVLEN'] = abs(int(variant['break1'])-int(variant['break2']))+1 # need to be an Integer, conversion is made afterwards
                rec = vcf_header.new_record(contig=str(variant['chrom1']),
                                                   start=int(variant['break1']),  # 1-based
                                                   stop=int(variant['break2']),  # 1-based
                                                   alleles=['N', '<DEL>'],
                                                   id='CallId'+str(i),
                                                   # ~ qual=float(variant['score']),
                                                   qual=30, # test value
                                                   filter=variant['filter'],
                                                   info=infos,
                                                   samples=samples)
                records.append(rec)               
                i += 1
            else:
                continue
    return records


def main(bedfile, output_file, genome_file):
    """
        Constructing a vcf file from scratch using the linkedSV
        bedpe inputfile thzta describes the variants
    """
    contigs = get_contigs(genome_file)
    tmp_vcf = vcf_from_scratch(output_file, contigs)
    vcf_in = VariantFile(tmp_vcf)

    ashkenazim_son = "HG002"

    vcf_in.header.info.add('END', number=1, type='Integer',
                           description="End position of the variant "
                                       "described in this record")
    vcf_in.header.info.add('SVLEN', number=1, type='Integer',
                           description="Length of the variant "
                                       "described in this record")
    vcf_in.header.info.add('SVTYPE', number=1, type='String',
                           description="Type of structural variant")
    vcf_in.header.info.add('SVMETHOD', number=1, type='String',
                           description="SV detection method")
    vcf_in.header.filters.add('FAIL', number=None, type=None,
                           description="Fail to pass filtering")                       
    vcf_in.header.formats.add('GT', number=1, type='String',
                              description="Genotype")
    vcf_in.header.add_sample(ashkenazim_son)

    records = reformat_bedpe2vcfrecords(bedfile, vcf_in.header)

    vcf_out = VariantFile(output_file, 'w', header=vcf_in.header)
    for rec in records:
        vcf_out.write(rec)

    os.remove(tmp_vcf)


def parse_arguments():
    parser = argparse.ArgumentParser(prog="naibr2vcf.py",
                                     description="Reformating NAIBR output")
    parser.add_argument("-i", "--input", required=True,
                        help="NAIBR bedpe file")
    parser.add_argument("-o", "--output", required=True,
                        help="NAIBR vcf output file")
    parser.add_argument("-r", "--reference", required=True,
                        help="The reference genome file previously "
                             "indexed with samtools faidx")
    args = parser.parse_args()
    return args.input, args.output, args.reference


if __name__ == "__main__":
    arguments = parse_arguments()
    main(*arguments)
