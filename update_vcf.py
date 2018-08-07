#!/usr/bin/env python

from __future__ import print_function, unicode_literals, absolute_import
from dirty_scrub import vcfcombiner
import argparse
import logging


def _parse_args():
    """
    Parse the input arguments.
    """
    parser = argparse.ArgumentParser(description='Flag VCF entries as clean')
    parser.add_argument('original_vcf', metavar='FILE', help='Input VCF file')
    parser.add_argument(
        'input_tsvs', metavar='FILES', nargs='*', help='TSVs to apply to input vcf.')
    parser.add_argument(
        "--output",
        type=str,
        help="output VCF file name. Defaults to prefix.scrubbed.vcf")

    return parser.parse_args()


if __name__ == '__main__':
    args = _parse_args()

    logging.basicConfig(filename='DIRTY_SCRUBBER.log', level=logging.DEBUG, filemode='w')

    # argparse
    if args.target_vcf is None:
        target_vcf = "{ori_vcf_prefix}.scrubbed.vcf".format(
            ori_vcf_prefix=args.original_vcf.rstrip(".vcf"))
    else:
        target_vcf = args.target_vcf

    # Combine flags
    with vcfcombiner.concat_read_only_file_stream(args.input_tsvs) as tsv_fhs:
        vcfcombiner.vcf_update(
            target_vcf, source_vcf=args.original_vcf, tsv_stream=tsv_fhs)  # TODO
