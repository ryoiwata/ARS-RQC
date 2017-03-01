#!/usr/bin/env python3
# rqcmain.py - The core module for ARS-RQC qulaity control workflow
# Adam Rivers 02/2017 USDA-ARS-GBRU
import argparse
import os
import subprocess
import logging
import tempfile
import json
import shutil


def build_vertebrate_db(cat, dog, mouse, human, datadir):
    """Builds a bbsplit.sh database for mapping reads to masked versions of
    the cat, dog, human and mouse genome. Returns the path of the database"""
    try:
        parameters = ['bbsplit.sh', 'build=1', 'k=14', 'usemodulo',
                      'ref_cat=' + os.path.abspath(cat),
                      'ref_dog=' + os.path.abspath(dog),
                      'ref_mouse=' + os.path.abspath(mouse),
                      'ref_mouse=' + os.path.abspath(human),
                      'path=' + os.path.join(datadir)]
        p0 = subprocess.run(parameters, check=True, stderr=subprocess.PIPE)
        return p0.stderr.decode('utf-8')
    except RuntimeError:
        print("Couldn't build DB of vertebrate contaminants using bbsplit")


def estimate_kmer_coverage(histogram, outdir):
    """estimates the proportion of the kmers at a depth of 3x, 5x and 10x \
    and extrapolates the coverage out beyond the current coverage using a \
    rational function approximation method in PreSeqR. Returns json object."""
    try:
        parameters = ['coverage_est.R',
                      '--plot', 'FALSE',
                      '--input', os.path.abspath(histogram),
                      '--coverage', 'c(3, 5, 10)']
        p6 = subprocess.run(parameters, check=True, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
        return p6.stdout.decode('utf-8')
    except RuntimeError:
        print("Could not estimate kmer coverage beyond current depth \
              using PreseqR")


class Fastq():

    @classmethod
    def parse_params(self):
        """converts parameters from a json file into a dictionary formatted
        for use by bbtools"""
        try:
            with open("data/parameters.json", 'r') as p:
                fastq_parameters = json.load(p)
            bbtoolsdict = {}
            for key in fastq_parameters["rqcfilter"]:
                bbtoolsdict[key] = []
                for key2 in fastq_parameters["rqcfilter"][key]:
                    if isinstance(fastq_parameters["rqcfilter"][key][key2],
                                  str):
                        bbtoolsdict[key].append(str(key2) + "=" +
                                                str(fastq_parameters
                                                ["rqcfilter"][key][key2]))
                    elif isinstance(fastq_parameters["rqcfilter"][key][key2],
                                    list):
                        bbtoolsdict[key].append(str(key2) + "=" +
                                                ",".join(fastq_parameters
                                                         ["rqcfilter"]
                                                         [key][key2]))
            return bbtoolsdict
        except RuntimeError:
                print("Could not load and parse the parameters.json file \
                      correctly")

    def __init__(self, path):
        self.abspath = os.path.abspath(path)
        self.filepath, self.filename = os.path.split(os.path.abspath(path))
        self.metadata = {}

    def __repr__(self):
        return 'Fastq Class object :' + self.filename

    def filter_contaminants(self, outdir):
        """Calls bbduk to perform adapter removal and create quality data"""
        try:
            bbtoolsdict = self.parse_params()
            parameters = ['bbduk.sh',
                          'in=' + self.abspath,
                          'out=' + os.path.join(outdir, 'clean1.fq.gz'),
                          # Write statistics about  contamininants detected.
                          'stats=' + os.path.join(outdir,
                                                  'scaffoldStats1.txt'),
                          # Base composition histogram by position
                          'bhist=' + os.path.join(outdir, 'bhist.txt'),
                          #  Quality histogram by position.
                          'qhist=' + os.path.join(outdir, 'qhist.txt'),
                          # Count of bases with each quality value.
                          'qchist=' + os.path.join(outdir, 'qchist.txt'),
                          # Histogram of average read quality.
                          'aqhist=' + os.path.join(outdir, 'aqhist.txt'),
                          # Quality histogram designed for box plots.
                          'bqhist=' + os.path.join(outdir, 'bqhist.txt'),
                          # Read GC content histogram.
                          'gchist=' + os.path.join(outdir, 'gchist.txt')
                          ]
            parameters.extend(bbtoolsdict['filter_contaminants'])
            p1 = subprocess.run(parameters, check=True, stderr=subprocess.PIPE)
            self.metadata['filter_contaminants'] = list(os.walk(outdir))
            return p1.stderr.decode('utf-8')
        except RuntimeError:
            print("could not perform contaminant filtering with bbduk")
            print(p1.stderr.decode('utf-8'))

    def trim_adaptors(self, outdir):
        """Calls bbduk to remove contaminant sequences"""
        try:
            bbtoolsdict = self.parse_params()
            parameters = ['bbduk.sh', 'in=' + self.abspath,
                          'out=' + os.path.join(outdir, 'clean2.fq.gz'),
                          'stats=' + os.path.join(outdir,
                                                  'scaffoldStats2.txt')]
            parameters.extend(bbtoolsdict['trim_adaptors'])
            p2 = subprocess.run(parameters, check=True, stderr=subprocess.PIPE)
            self.metadata['trim_adaptors'] = list(os.walk(outdir))
            return p2.stderr.decode('utf-8')
        except RuntimeError:
            print("could not perform adaptor removal with bbduk")

    def merge_reads(self, outdir):
        """merge reads to generate insert size histogram and \
        error correct, retuns unmerged reads"""
        try:
            bbtoolsdict = self.parse_params()
            parameters = ['bbmerge.sh', 'in1=' + self.abspath,
                          'ihist=' + os.path.join(outdir,
                                                  'merge_histogram.txt'),
                          'outc=' + os.path.join(outdir, 'cardinality.txt'),
                          'outm=' + os.path.join(outdir, 'merged.fq.gz'),
                          'outu=' + os.path.join(outdir, 'unmerged.fq.gz')]
            parameters.extend(bbtoolsdict['merge_reads'])
            p3 = subprocess.run(parameters, stderr=subprocess.PIPE)
            self.metadata['merge_reads'] = list(os.walk(outdir))
            return p3.stderr.decode('utf-8')
        except RuntimeError:
            print("could not perform read merging with bbmerge")

    def remove_vertebrate_contaminants(self, outdir):
        """maps reads to repeat-masked human, dog, cat and mouse genomes
        to remove contaminants"""
        try:
            bbtoolsdict = self.parse_params()
            parameters = ['bbsplit.sh',
                          'in=' + self.abspath,
                          'outu=' + os.path.join(outdir, 'novert.fq.gz')]
            parameters.extend(bbtoolsdict['remove_vertebrate_contaminants'])
            p4 = suborocess.run(parameters, check=True, stderr=subprocess.PIPE)
            self.metadata['remove_vertebrate_contaminants'] = list(
                          os.walk(outdir))
            return p4.stderr.decode('utf-8')
        except RuntimeError:
            print("Could not perform vertebrate conaminant removal with bbmap")

    def sortbyname(self):
        """Sorts a fastq file by read names, outputs uncompressed fastq"""
        try:
            temp_ordered_dir = tempfile.mkdtemp()
            root, base = os.path.split(self.abspath)
            if base.endswith(".gz"):
                uc = base.rsplit(".", 1)[0]
            else:
                uc = base
            tfile = os.path.join(temp_ordered_dir, 'sorted.fq')
            parameters = ['sortbyname.sh', 'in=' + self.abspath,
                          'out=' + tfile]
            p5 = subprocess.run(parameters, check=True, stderr=subprocess.PIPE)
            shutil.move(tfile, os.path.join(root, uc))
            return p5.stderr.decode('utf-8')
        except RuntimeError:
            print("could not reorder fastq file by name")
            print(p5.stderr.decode('utf-8'))
        finally:
            shutil.rmtree(temp_ordered_dir)

    def calculate_kmer_histogram(self, outdir):
        """calcualtes kmer histogram from a fastq file using BBtools \
        kmercountexact.sh and if that fails, approximates it with khist.sh"""
        try:
            parameters = ['kmercountexact.sh',
                          'in=' + self.abspath,
                          'hist=' + os.path.join(outdir, 'kmerhist.txt')]
            p5a = subprocess.run(parameters, check=True,
                                 stderr=subprocess.PIPE)
            self.metadata['calculate_kmer_histogram'] = list(os.walk(outdir))
            return p5a.stderr.decode('utf-8')
        except RuntimeError:
            print("Could not calculate the kmer histogram with \
                  kmerexactcount.sh. Attempting to estimate it with khist.sh")
            try:
                parameters = ['khist.sh',
                              'in=' + self.abspath,
                              'histcol=2',
                              'hist=' + os.path.join(outdir, 'kmerhist.txt')]
                p5b = subprocess.run(parameters, check=True,
                                     stderr=subprocess.PIPE)
                self.metadata['calculate_kmer_histogram'] = list(
                              os.walk(outdir))
            except RuntimeError:
                return p5b.stderr.decode('utf-8')

    def clumpify(self, outdir):
        """Reorders reads or read pairs in a fastq file by shared kmers. \
        This reduces the size of compressed files by about 30% \
        and speeds up kmer-based analyses like de Bruijn assembly by \
        increasing the use of CPU cashe"""
        try:
            bbtoolsdict = self.parse_params()
            parameters = ['clumpify.sh',
                          'in=' + self.abspath,
                          'out=' + os.path.join(outdir, 'clumped.fq.gz')]
            parameters.extend(bbtoolsdict['clumpify'])
            p6 = subprocess.run(parameters, check=True,
                                stderr=subprocess.PIPE)
            self.metadata['clumpify'] = list(os.walk(outdir))
            return p6.stderr.decode('utf-8')
        except RuntimeError:
            print("Could not reorder and error correct the data with \
                  clumpify.sh")
            return p6.stderr.decode('utf-8')
