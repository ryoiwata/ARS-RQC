import unittest
import os
import tempfile
import filecmp
import sys
import rqcfilter
import shutil

class TestfastqMethods(unittest.TestCase):

    def test_filter_contaminants(self):
        testdir = tempfile.TemporaryDirectory()  # Create temp dir
        obj = rqcfilter.Fastq('data/reads.fq.gz')  # create fastq object
        obj.filter_contaminants(os.path.join(testdir.name))  # run method
        fq1 = os.path.join(testdir.name, 'clean1.fq.gz')  # test output datadir
        fqo1 = rqcfilter.Fastq(fq1)  # create fastq object from output fastq
        fqo1.sortbyname()  # sort by read name (read order from bbduk varies)
        fq2 = 'tests/outputs/filter_contaminants/clean1.fq.gz'  # valid. data
        fqo2 = rqcfilter.Fastq(fq2)  # create object with validation data
        fqo2.sortbyname()  # sort validation data by name
        shutil.copytree(testdir.name, 'fc_dir')
        self.assertTrue(filecmp.cmp(fq1, fq2))  # test the outputs are the same

    def test_trim_adaptors(self):
        testdir2 = tempfile.TemporaryDirectory()  # Create temp dir
        obj = rqcfilter.Fastq('data/reads.fq.gz')  # create fastq object
        obj.trim_adaptors(testdir2.name)  # run method
        fq1 = os.path.join(testdir2.name, 'clean2.fq.gz')  # test output ddir
        fqo1 = rqcfilter.Fastq(fq1)  # create fastq object from output fastq
        fqo1.sortbyname()  # sort by read name (read order from bbduk varies)
        fq2 = 'tests/outputs/trim_adaptors/clean2.fq.gz'  # valid. data
        fqo2 = rqcfilter.Fastq(fq2)  # create object with validation data
        fqo2.sortbyname()  # sort validation data by name
        shutil.copytree(testdir2.name, 'ta_dir')
        self.assertTrue(filecmp.cmp(fq1, fq2))  # test the outputs are the same

    def test_merge_reads(self):
        testdir3 = tempfile.TemporaryDirectory()  # Create temp dir
        print(testdir3.name)
        obj = rqcfilter.Fastq('data/reads.fq.gz')  # create fastq object
        obj.merge_reads(testdir3.name)  # run method
        # sort merged reads from test
        fqm1 = os.path.join(testdir3.name, 'merged.fq.gz')
        fqmo1 = rqcfilter.Fastq(fqm1)
        fqmo1.sortbyname()
        #sort unmerged reads from test
        fqu1 = os.path.join(testdir3.name, 'unmerged.fq.gz')
        fquo1 = rqcfilter.Fastq(fqu1)
        fquo1.sortbyname()
        # sort merged reads from validation data
        fqm2 = 'tests/outputs/trim_adaptors/merged.fq.gz'  # valid. data
        fqmo2 = rqcfilter.Fastq(fqm2)  # create object with validation data
        fqmo2.sortbyname()  # sort validation data by name
        # sort unmerged reads from validation data
        fqu2 = 'tests/outputs/trim_adaptors/unmerged.fq.gz'  # valid. data
        fquo2 = rqcfilter.Fastq(fqu2)  # create object with validation data
        fquo2.sortbyname()  # sort validation data by name
        # compare files
        shutil.copytree(testdir.name, 'mr_dir')
        self.assertTrue(filecmp.cmp(fqm1, fqm2) and filecmp.cmp(fqu1, fqu2))


if __name__ == '__main__':
    unittest.main()
