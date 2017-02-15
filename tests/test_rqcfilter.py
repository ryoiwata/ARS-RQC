import unittest
import rqcfilter
import os
import tempfile


class TestfastqMethods(unittest.TestCase):

    def test_trim_adaptors(self):
        with open(sample_trim_adaptors.txt, 'r') as sta:
            sta_string = sta.readlines()
        testdir = tempfile.tempdir()
        rqc_string = rqcfilter.Fastq.trim_adaptors(os.path.abspath
                                                   ('../data/reads100.fq.gz'),
                                                   testdir)
        self.assertEqual(sta_string, rqc_string)


if __name__ == '__main__':
    unittest.main()
