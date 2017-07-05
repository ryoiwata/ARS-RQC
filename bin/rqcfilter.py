#!/usr/bin/env python3
# rqcfilter.py - A sequence quality control and metadata collection workflow
# Adam Rivers 02/2017 USDA-ARS-GBRU
import argparse
import os
import subprocess
import shutil
import logging
import tempfile
import json
import shutil
import sys
import pandas as pd
import numpy as np
from ars_rqc import rqcmain
from ars_rqc import rqcparser
from ars_rqc.definitions import ROOT_DIR


# Utility functions

# this extends the json class to handle numpy types in dictionaries
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(MyEncoder, self).default(obj)

def convert_keys_to_string(dictionary):
    """Recursively converts dictionary keys to strings."""
    if not isinstance(dictionary, dict):
        return dictionary
    return dict((str(k), convert_keys_to_string(v))
                for k, v in dictionary.items())


def write_metadata(indir, outfile):
    """Writes dictionary to json file"""
    try:
        datadict = rqcparser.parse_dir(indir)
    except RuntimeError:
        print("Could not parse bbtools output file(s)")
    try:
        with open(outfile, 'w') as fp:
            json.dump(datadict, fp, cls=NumpyEncoder)
    except IOError:
        print("Could not write json metadata file")


def create_clean_name(fastq):
    """Parses the fastq input name, inserting 'arsrqc' into the name"""
    try:
        fastqnamelist = os.path.basename(fastq).split(".")
        if fastqnamelist[-1] is 'gz' and fastqnamelist[-2] in ('fastq', 'fq'):
            rqcnamelist = fastqnamelist.insert(-2, 'arsrqc')
            rqcname = '.'.join(rqcnamelist)
        elif fastqnamelist[-1] in ('fastq', 'fq'):
            rqcnamelist = fastqnamelist.insert(-1, 'arsrqc')
            rqcnamelist.append('gz')
            rqcname = '.'.join(rqcnamelist)
            return rqcname
    except RuntimeError:
        print("Could not parse the name of the input fastq file. Please \
        make sure it ends in .fq, .fastq, .fq.gz or .fastq.gz ")


def mk_temp_dir(tempdir, suffix):
    """takes a root directory and a suffix and creates that \
    directory, logging progress and returning the path"""
    fulltemp = os.path.join(tempdir, suffix)
    logging.info('Creating temporary directory: {}'.format(fulltemp))
    try:
        os.mkdir(fulltemp)
    except IOError:
        logging.error('could not create temporary \
                      directory: {}'.format(fulltemp))
    return fulltemp


def myparser():
    parser = argparse.ArgumentParser(description='rqcfilter.py - \
                                     A sequence quality control and metadata \
                                     collection workflow.')

    parser.add_argument('--fastq', '-f', type=str, required=True,
                        help='A .fastq, .fq, .fastq.gz or .fq.gz file.')

    parser.add_argument('--output', '-o', type=str, default='rqcout',
                        help='the output directory')
    parser.add_argument('--overwrite', '-w', action='store_true', default=False,
                        help='a flag to overwrite the output directory')

    parser.add_argument('--removevertebrates', '-r', action='store_true',
                        default=False,
                        help='A flag to specify if the reads should be mapped  \
                        against human, cat, dog and mouse genomes to find \
                        contaminants. Default is false.')

    parser.add_argument('--paired', '-p', action='store_true', default=False,
                        help='A flag to specify if the fastq file is \
                        paired and interleaved. Default is false.')

    parser.add_argument('--keepmergeresults', '-m', action='store_true',
                        default=False, help='A flag to specify whether to keep \
                        a fastq file with merged reads and a fastq file with \
                        umerged reads. Defalt is false')
    parser.add_argument('--keepfullresults', '-k', action='store_true',
                        default=False, help='A flag to specify whether to keep \
                        all intermediate files or just the summary log, \
                        sequence and metadata files')
    args = parser.parse_args()
    return args


# Main loop
def main():

    args = myparser()  # load command line options

    # Create the output directory
    if args.overwrite:  # if overwrite is true:
        if os.path.exists(args.output):
            shutil.rmtree(args.output)
            os.makedirs(args.output)
        else:
            os.makedirs(args.output)
    else:  # if overwrite is false:
        if os.path.exists(args.output):
            raise IOError('The ouput directory already exists, check the \
                          output flag')
        else:
            os.makedirs(args.output)

    # Create temporary directory
    rqctempdir = tempfile.mkdtemp()

    # Set up logging
    logging.basicConfig(filename=os.path.join(args.output, 'rqc.log'),
                        level=logging.INFO,
                        format='%(asctime)s %(message)s')
    logging.info('Starting USDA ARS GBRU rolling quality control workflow.')

    # Assign globals
    abs_fastq = os.path.abspath(args.fastq)

    # Remove contaminants
    tmp_fc = mk_temp_dir(rqctempdir, 'filter_contaminants')  # make temp. dir.
    logging.info('Starting contaminant removal')
    data1 = rqcmain.Fastq(path=abs_fastq)  # create Fastq object
    decondata = data1.filter_contaminants(outdir=tmp_fc)  # decontaminate
    logging.info(decondata)  # record bbduk output in log

    # Trim adaptors
    tmp_ta = mk_temp_dir(rqctempdir, 'trim_adaptors')  # make temp. dir.
    logging.info('Starting adaptor trimming')
    data2 = rqcmain.Fastq(os.path.join(tmp_fc, 'clean1.fq.gz'))  # create fastq
    trimdata = data2.trim_adaptors(outdir=tmp_ta)  # trim adaptors
    logging.info(trimdata)  # record output of second bbduk run

    # Remove vertebrate contaminants
    if args.removevertebrates:
        tmp_rvc = mk_temp_dir(rqctempdir, 'remove_vertebrate_contaminants')
        logging.info('Removing dog, cat, mouse and human reads')
        # Create new Fastq object
        data3 = rqcmain.Fastq(os.path.join(tmp_ta, 'clean2.fq.gz'))
        if not os.path.isdir(data/dogcatmousehuman/ref):
            rqcmain.build_vertebrate_db(cat=data/cat.fa.gz,
                                        dog=data/dog.fa.gz,
                                        human=data/hg19.fa.gz,
                                        mouse=data/mouse.fa.gz,
                                        datadir=data/dogcatmousehuman)
        rvcdata = data3.remove_vertebrate_contaminants(outdir=tmp_rvc)
        logging.info(rvcdata)

    # Clumpify data (Order reads by overlapping kmers to increase \
    # compression and processing speed)
    tmp_cy = mk_temp_dir(rqctempdir, 'clumpify')  # make temp dir
    logging.info('Starting contaminant removal')
    if args.removevertebrates:  # Select input file
        infile = os.path.join(tmp_rvc, 'novert.fq.gz')
    else:
        infile = os.path.join(tmp_ta, 'clean2.fq.gz')
    data4 = rqcmain.Fastq(path=infile)  # create Fastq object
    clumpdata = data1.clumpify(outdir=tmp_cy)  # clumpify
    logging.info(clumpdata)

    # Merge files
    if args.paired:
        tmp_mr = mk_temp_dir(rqctempdir, 'merge_reads')  # make temp. dir.
        logging.info('Merging read pairs')
        data5 = rqcmain.Fastq(os.path.join(tmp_cy, 'clumped.fq.gz'))
        mergedata = data5.merge_reads(outdir=tmp_mr)
        logging.info(mergedata)

    # Calculate kmer histogram
    tmp_kh = mk_temp_dir(rqctempdir, 'calculate_kmer_histogram')
    logging.info("calculating Kmer Histogram")
    khdata = data5.calculate_kmer_histogram(outdir=tmp_kh)
    logging.info(khdata)

    # TODO
    # Run PreseqR once the interface is setup and the R script has been fixed

    # copy all files from the tempdir to the output dir
    if args.keepfullresults:
        # copy all files from the tempdir to the output dir
        try:
            logging.info("Copying files from temporary directory to ouput \
                         directory")
            shutil.copytree(rqctempdir, os.path.join(args.output, "output"))
            logging.info("Parsing the metadata and writing it to a json file")
            write_metadata(indir=rqctempdir,
                           outfile=os.path.join(args.output, 'metadata.json'))
            shutil.rmtree(rqctempdir)
        except RuntimeError:
            print("could not copy the temproary directory to the ")
    # Create name for clean file
    else:
        try:
            cleanname = create_clean_name(args.fastq)
            rqcloc = os.path.join(args.output, cleanname)
            logging.info("Copying RQC processed fastq to {}".format(rqcloc))
            shutil.copy2(os.path.join(tmp_cy, 'clumped.fq.gz'), rqcloc)
            logging.info("Parsing the metadata and writing it to a json file")
            write_metadata(indir=rqctempdir,
                           outfile=os.path.join(args.output, 'metadata.json'))
            if args.keepmergeresults:
                cm1 = cleanname.split('.')
                cmm = fastqnamelist.insert(-2, 'merged')
                cmms = '.'.join(cmm)
                cmu = fastqnamelist.insert(-2, 'unmerged')
                cmus = '.'.join(cmus)
                shutil.copy2(os.path.join(tmp_mr, 'merged.fq.gz'),
                             os.path.join(args.output, cmms))
                shutil.copy2(os.path.join(tmp_mr, 'unmerged.fq.gz'),
                             os.path.join(args.output, cmus))
        except RuntimeError:
            print("Could not move all files to the output directory.")

    logging.info("Completed RQC run")


if __name__ == '__main__':
    main()
