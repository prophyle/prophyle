#! /usr/bin/env python3

import sys, unittest
from itertools import zip_longest

class exk_test(unittest.TestCase):

    def test_results(self):
        with open("expected_results.txt", 'r') as exp_results, open("results.txt", 'r') as results:
            for e_line, r_line in zip_longest(exp_results, results):
                self.assertEqual(e_line, r_line)

if __name__ == '__main__':
    unittest.main()
