#!/usr/bin/env python
from __future__ import unicode_literals, print_function, division, absolute_import
from dirty_scrub.vcfcombiner import concat_read_only_file_stream, load_vcf_to_dict, vcf_update
from dirty_scrub.declarations import YomoDict

import unittest
import filecmp
import json


class TestVCFCombiner(unittest.TestCase):
    def test_concat_file_stream(self):
        test_file_1 = 'test_data/samp1.tsv'
        test_file_2 = 'test_data/samp2.tsv'
        expected_word = '10\t869444\n9921\t2103219\n21\t34\n22231_2\t6402439\n3332_GG\t100000\n'
        actual_word = ''
        with concat_read_only_file_stream(test_file_1, test_file_2) as f:
            for line in f:
                actual_word += line
        self.assertEqual(expected_word, actual_word)

    def test_load_vcf(self):
        with open('test_data/H002.lumpy_noalts.json', 'r') as f:
            expected_dict = json.load(f)

        actual_dict = load_vcf_to_dict('test_data/H002.lumpy_noalts.vcf')

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


class DirtyScrubberIntegration(unittest.TestCase):
    #@unittest.skip("Not ready yet")
    def test_vcf_flag(self):
        # Setup
        source_vcf = 'test_data/H002.lumpy_noalts.vcf'
        expected_flagged_vcf = 'test_data/H002.lumpy_noalts.TEST-FLAGGED.vcf'
        actual_flagged_vcf = 'DELETEME.vcf'
        test_file_1 = 'test_data/samp1.tsv'
        test_file_2 = 'test_data/samp2.tsv'

        # do it
        with concat_read_only_file_stream(test_file_1, test_file_2) as tsv_stream:
            vcf_update(
                target_vcf=actual_flagged_vcf, source_vcf=source_vcf, tsv_stream=tsv_stream)

        self.assertTrue(
            filecmp.cmp(expected_flagged_vcf, actual_flagged_vcf))


if __name__ == '__main__':
    unittest.main(verbosity=2)
