#!/usr/env/python3
# test_rqcmain.py - a testing module for rqcmain.py
# Adam Rivers 02/2017 USDA-ARS-GBRU

import unittest
import os
import tempfile
import filecmp
import sys
from ars_rqc import rqcmain
import shutil

class TestfastqMethods(unittest.TestCase):

    def test_filter_contaminants(self):
        try:
            testdir = tempfile.mkdtemp()
            obj = rqcmain.Fastq('data/reads.fq.gz')  # create fastq object
            obj.filter_contaminants(testdir)  # run method
            fq1 = os.path.join(testdir, 'clean1.fq.gz')  # test output datadir
            fqo1 = rqcmain.Fastq(fq1)  # create fastq object from output fastq
            fqo1.sortbyname()  # sort by read name (read order from bbduk varies)
            fq2 = 'tests/outputs/filter_contaminants/clean1.fq'
            fqo2 = rqcmain.Fastq(fq2)  # create object with validation data
            fqo2.sortbyname()  # sort validation data by name
            #shutil.copytree(testdir, 'fc_dir')
            self.assertTrue(filecmp.cmp(os.path.join(testdir, 'clean1.fq'), fq2))  # test the outputs are the same
        except RuntimeError:
            print("could not compare results of filter contaminants")
        finally:
            shutil.rmtree(testdir)


    def test_trim_adaptors(self):
        try:
            testdir2 = tempfile.mkdtemp()
            obj = rqcmain.Fastq('data/reads.fq.gz')  # create fastq object
            obj.trim_adaptors(testdir2)  # run method
            fq1 = os.path.join(testdir2, 'clean2.fq.gz')  # test output ddir
            fqo1 = rqcmain.Fastq(fq1)  # create fastq object from output fastq
            fqo1.sortbyname()  # sort by read name (read order from bbduk varies)
            fq2 = 'tests/outputs/trim_adaptors/clean2.fq'  # valid. data
            fqo2 = rqcmain.Fastq(fq2)  # create object with validation data
            fqo2.sortbyname()  # sort validation data by name
            self.assertTrue(filecmp.cmp(os.path.join(testdir2, "clean2.fq"), fq2))  # test the outputs are the same
        except RuntimeError:
            print("could not compare results of filter contaminants")
        finally:
            shutil.rmtree(testdir2)

    def test_merge_reads(self):
        try:
            testdir3 = tempfile.mkdtemp()
            obj = rqcmain.Fastq('data/reads.fq.gz')  # create fastq object
            obj.merge_reads(testdir3)  # run method
            # sort merged reads from test
            fqm1 = os.path.join(testdir3, 'merged.fq.gz')
            fqmo1 = rqcmain.Fastq(fqm1)
            fqmo1.sortbyname()
            # sort unmerged reads from test
            fqu1 = os.path.join(testdir3, 'unmerged.fq.gz')
            fquo1 = rqcmain.Fastq(fqu1)
            fquo1.sortbyname()
            # sort merged reads from validation data
            fqm2 = 'tests/outputs/merge_reads/merged.fq'  # valid. data
            fqmo2 = rqcmain.Fastq(fqm2)  # create object with validation data
            fqmo2.sortbyname()  # sort validation data by name
            # sort unmerged reads from validation data
            fqu2 = 'tests/outputs/merge_reads/unmerged.fq'  # valid. data
            fquo2 = rqcmain.Fastq(fqu2)  # create object with validation data
            fquo2.sortbyname()  # sort validation data by name
            # compare files
            self.assertTrue(filecmp.cmp(os.path.join(testdir3, "merged.fq"),
                            fqm2) and filecmp.cmp(os.path.join(testdir3,
                            "unmerged.fq"), fqu2))
        except RuntimeError:
            print("could not compare results of filter contaminants")
        finally:
            shutil.rmtree(testdir3)


if __name__ == '__main__':
    unittest.main()
