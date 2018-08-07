#!/usr/bin/env python
from __future__ import unicode_literals, print_function, division, absolute_import
from dirty_scrub.vcfcombiner import concat_read_only_file_stream

import unittest


class TestQualiMap2(unittest.TestCase):
    def test_concat_file_stream(self):
        test_file_1 = 'test_data/samp1'
        test_file_2 = 'test_data/samp2'
        expected_word = 'Hello\nWorld\nBye\nWorld\n'
        actual_word = ''
        with concat_read_only_file_stream(test_file_1, test_file_2) as f:
            for line in f:
                actual_word += line
        self.assertEquals(expected_word, actual_word)


if __name__ == '__main__':
    unittest.main()
