#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
  Converting pLinkedSV bedpe to vcf format

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
    linkedsv_variants = BedTool(bedpefile)
    # Bedpe fields : https://github.com/WGLab/LinkedSV
    keys = ['chrom1', 'start1', 'end1',
            'chrom2', 'start2', 'end2',
            'sv_type', 'sv_id', 'sv_length',
            'qual', 'filter', 'info']

    records = []
    infos = defaultdict()
    for row in linkedsv_variants:
        variant = dict(zip(keys, row))
        infos.clear()
        infos['SVTYPE'] = variant['sv_type']
        for info in variant['info'].split(","):
            key, value = info.split('=')
            infos[key] = value
        samples = [{'GT': (None, None)}]
        rec = vcf_header.new_record(contig=str(variant['chrom1']),
                                           start=int(variant['end1']),  # 1-based
                                           stop=int(variant['end2']),  # 1-based
                                           alleles=['N', '<DEL>'],
                                           id=variant['sv_id'],
                                           qual=int(variant['qual']),
                                           filter=variant['filter'],
                                           info=infos,
                                           samples=samples)
        records.append(rec)
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
    vcf_in.header.info.add('SVTYPE', number=1, type='String',
                           description="Type of structural variant")
    vcf_in.header.info.add('SVMETHOD', number=1, type='String',
                           description="SV detection method")
    vcf_in.header.formats.add('GT', number=1, type='String',
                              description="Genotype")
    vcf_in.header.add_sample(ashkenazim_son)

    records = reformat_bedpe2vcfrecords(bedfile, vcf_in.header)

    vcf_out = VariantFile(output_file, 'w', header=vcf_in.header)
    for rec in records:
        vcf_out.write(rec)

    os.remove(tmp_vcf)


def parse_arguments():
    parser = argparse.ArgumentParser(prog="filter.py",
                                     description="Reformating linkedSV output")
    parser.add_argument("-i", "--input", required=True,
                        help="LinkedSV bedpe file")
    parser.add_argument("-o", "--output", required=True,
                        help="LinkedSV vcf output file")
    parser.add_argument("-r", "--reference", required=True,
                        help="The reference genome file previously "
                             "indexed with samtools faidx")
    args = parser.parse_args()
    return args.input, args.output, args.reference


if __name__ == "__main__":
    arguments = parse_arguments()
    main(*arguments)
