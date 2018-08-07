#!/usr/bin/env python

from __future__ import print_function, unicode_literals, absolute_import
from itertools import chain
from contextlib import contextmanager
import logging
import csv
import os

from .declarations import YomoDict


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
        filter_codes = [elem.strip() for elem in vcf_row[6].split(';')]

        if variant_id == '.':
            logging.debug('Encountered "." on line {line_num}'.format(
                line_num=line_num))

        try:
            vcf_obj[variant_id] = {'POS': int(variant_start_pos), 'filter_codes': filter_codes}
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
    return var_id in vcf_obj and vcf_obj[var_id]['POS'] == var_pos


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

    TODO:
        * Refactor to:
            * open source_vcf only once, load VCF into relevant data structure
            * Update vcf obj with tsv stream (in parallizable manner)
            * Write to VCF, retaining all fields and comments
              * Add a FILTER comment describt the DIRT filter code
    """
    def flag_vcf_entry(filters):
        'DIRT' not in filters and filters.append('DIRT')
        logging.debug('New filters list\n' + str(filters))

    ori_vcf_obj = load_vcf_to_dict(source_vcf)
    tsv_reader = csv.reader(tsv_stream, delimiter='\t')

    for tsv_row in tsv_reader:
        logging.debug("TSV ROW: " + str(tsv_row))
        var_id, var_pos = tsv_row[0], int(tsv_row[1])
        if should_flag_variant(ori_vcf_obj, var_id=var_id, var_pos=var_pos):
            logging.info("FLAGGING: {var_id}\t{var_pos}".format(var_id=var_id, var_pos=var_pos))
            flag_vcf_entry(ori_vcf_obj[var_id]['filter_codes'])

    with open(target_vcf, 'w') as new_vcf, open(source_vcf, 'r') as old_vcf:
        tsv_reader = csv.reader(old_vcf, delimiter='\t')
        tsv_writer = csv.writer(new_vcf, delimiter='\t')

        for row in tsv_reader:
            var_id = row[2]
            if len(row) > 8 and var_id in ori_vcf_obj:
                filter_codes = ';'.join(ori_vcf_obj[var_id]['filter_codes'])
                logging.debug('Filter codes: ' + str(filter_codes))
                row[6] = filter_codes
                logging.debug('ALTERED ROW:\n' + str('\t'.join(row)))

            tsv_writer.writerow(row)
