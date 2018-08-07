#!/usr/bin/env python

from __future__ import print_function, unicode_literals, absolute_import
from itertools import chain
from contextlib import contextmanager
import logging
import csv

from .declarations import YomoDict
import pprint


FIXED_VCF_COLUMNS = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]


def vcf_entry_match(vcf_id, start, stop, target_vcf):
    """Returns whether the requested update entry matches the target VCF entry"""
    pass


@contextmanager
def concat_read_only_file_stream(*file_names):
    file_handlers = []
    try:
        for fn in file_names:
            file_handler = open(fn, 'r')
            file_handlers.append(file_handler)

        yield chain(*file_handlers)
    finally:
        for fh in file_handlers:
            fh.close()


def load_vcf_to_dict(vcf_file):
    """Loads VCF entries to a mapping ID -> POS

    expects row to be at least size 8 (VCF spec)
    """
    def _update_yomo_dict(vcf_row, line_num):
        logging.debug('Updating line: ' + str(line_num))
        variant_id = vcf_row[2]
        variant_start_pos = vcf_row[1]

        if variant_id == '.':
            logging.debug('Encountered "." on line {line_num}'.format(
                line_num=line_num))

        try:
            vcf_obj[variant_id] = int(variant_start_pos)
        except KeyError:
            logging.info('{var_id} id occured multiple times'.format)

    vcf_obj = YomoDict()
    with open(vcf_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for assumed_line, row in enumerate(reader):
            if len(row) < 8:
                logging.debug("Skipping line: " + str(assumed_line))
                continue
            _update_yomo_dict(vcf_row=row, line_num=assumed_line)

    logging.debug("LOADED VCF {vcf_file}:\n{vcf_obj}".format(
        vcf_file=vcf_file, vcf_obj=vcf_obj))

    return vcf_obj


def should_flag_variant(vcf_obj, var_id, var_pos):
    return var_id in vcf_obj and vcf_obj[var_id] == var_pos


def vcf_update(target_vcf, source_vcf, tsv_stream):
    """Top level functions for updating the target vcf FILTER column

    Strategy:
        1. Load VCF into memory
            * TODO benchmark this versus file handles
        2. Create data structure for quickly querying ID POS  # Just a mapping ID -> POS
            * TODO maybe pandas or some library may be better
        3. Stream in input TSVs
        4. Update for each row in input stream
        5. Write updated VCF to disk (target file)
    """
    ori_vcf_obj = load_vcf_to_dict(source_vcf)
    tsv_reader = csv.reader(tsv_stream, delimiter='\t')

    with open(target_vcf, 'w') as vcf:
        for tsv_row in tsv_reader:
            logging.debug("TSV ROW: " + str(tsv_row))
            var_id, var_pos = tsv_row[0], int(tsv_row[1])
            if should_flag_variant(ori_vcf_obj, var_id=var_id, var_pos=var_pos):
                logging.info("FLAGGING: {var_id}\t{var_pos}".format(var_id=var_id, var_pos=var_pos))
                vcf.write('FLAGGED {var_id}\n'.format(var_id=var_id))
