#!/usr/bin/env python
from __future__ import unicode_literals, print_function, division, absolute_import
from dirty_scrub.vcfcombiner import concat_read_only_file_stream, load_vcf_to_dict
from dirty_scrub.declarations import YomoDict

import unittest
import json


class TestVCFCombiner(unittest.TestCase):
    def test_concat_file_stream(self):
        test_file_1 = 'test_data/samp1'
        test_file_2 = 'test_data/samp2'
        expected_word = 'Hello\nWorld\nBye\nWorld\n'
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


if __name__ == '__main__':
    unittest.main(verbosity=2)
