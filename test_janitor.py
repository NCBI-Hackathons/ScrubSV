#!/usr/bin/env python
from __future__ import unicode_literals, print_function, division, absolute_import
from janitor.vcfcombiner import concat_read_only_file_stream, load_vcf_to_dict, vcf_update
from janitor.declarations import YomoDict

import unittest
import logging
import json
import os

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
TEST_DATA_FLD = os.path.join(SCRIPT_DIR, 'test_data')


def compare_file_content(file1, file2):
    with open(file1, 'r') as f1, open(file2, 'r') as f2:
        for line1, line2 in zip(f1, f2):
            if line1 != line2:
                return False
    return True


class TestJanitorTools(unittest.TestCase):
    def test_concat_file_stream(self):
        test_file_1 = os.path.join(TEST_DATA_FLD, 'samp1.tsv')
        test_file_2 = os.path.join(TEST_DATA_FLD, 'samp2.tsv')
        expected_word = '10\t869444\n9921\t2103219\n21\t34\n22231_2\t6402439\n3332_GG\t100000\n'
        actual_word = ''
        with concat_read_only_file_stream(test_file_1, test_file_2) as f:
            for line in f:
                actual_word += line
        self.assertEqual(expected_word, actual_word)

    def test_load_vcf(self):
        variant_json = os.path.join(TEST_DATA_FLD, 'H002.lumpy_noalts.json')
        variant_vcf = os.path.join(TEST_DATA_FLD, 'H002.lumpy_noalts.vcf')
        with open(variant_json, 'r') as f:
            expected_dict = json.load(f)

        actual_dict = load_vcf_to_dict(variant_vcf)

        self.assertEqual(expected_dict, actual_dict)


class TestYOMODict(unittest.TestCase):
    def setUp(self):
        self.test_dict = YomoDict()

    def test_only_add_key_once(self):
        self.test_dict['blah'] = 5

        try:
            self.test_dict['blah'] = 10
        except KeyError:
            return
        else:
            raise self.fail("Succeeded updating key")


class JanitorIntegration(unittest.TestCase):
    def test_vcf_flag(self):
        # log file
        logging.basicConfig(filename='TEST_janitor.log', level=logging.DEBUG, filemode='w')

        # Setup
        source_vcf = os.path.join(TEST_DATA_FLD, 'H002.lumpy_noalts.vcf')
        expected_flagged_vcf = os.path.join(TEST_DATA_FLD, 'H002.lumpy_noalts.TEST-FLAGGED.vcf')
        actual_flagged_vcf = os.path.join(TEST_DATA_FLD, 'DELETEME.vcf')
        test_file_1 = os.path.join(TEST_DATA_FLD, 'samp1.tsv')
        test_file_2 = os.path.join(TEST_DATA_FLD, 'samp2.tsv')

        # do it
        with concat_read_only_file_stream(test_file_1, test_file_2) as tsv_stream:
            vcf_update(
                target_vcf=actual_flagged_vcf, source_vcf=source_vcf,
                tsv_stream=tsv_stream)

        self.assertTrue(
            compare_file_content(expected_flagged_vcf, actual_flagged_vcf))


if __name__ == '__main__':
    unittest.main(verbosity=2)
