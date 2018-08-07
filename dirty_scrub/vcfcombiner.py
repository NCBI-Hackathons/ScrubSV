#!/usr/bin/env python

from __future__ import print_function, unicode_literals, absolute_import
from itertools import chain
from contextlib import contextmanager
import csv


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
    vcf_obj = dict()
    with open(vcf_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if len(row) < 8:
                continue
            variant_id = row[2]
            variant_start_pos = row[1]
            vcf_obj[variant_id] = int(variant_start_pos)

    return vcf_obj


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
    pass
