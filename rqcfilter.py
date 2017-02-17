#!/usr/env/python3
# rqcfilter.py - A sequence quality control and metadata collection workflow
# Adam Rivers 02/2017 USDA-ARS-GBRU
import argparse
import os
import subprocess
import logging
import tempfile
import json
import shutil
import rqcmain

def myparser():
    parser = argparse.ArgumentParser(description='rqcfilter.py - \
                                     A sequence quality control and metadata \
                                     collection workflow.')

    parser.add_argument('--fastq', '-f', type=str, required=True,
                        help='A .fastq, .fq, .fastq.gz or .fq.gz file.')

    parser.add_argument('--metadata.yaml', '-m', type=argparse.FileType('r'),
                        help='A yaml formated metadata file. \
                        See documentation.')

    parser.add_argument('--qcdata', '-qc', type=str, default='qcdata.json',
                        help='The name of the file containing QC data in \
                        json format')

    parser.add_argument('--output', '-o', type=str, default='rqcout',
                        help='the output directory')

    parser.add_argument('--removevertebrates', '-r', action='store_true',
                        default=False,
                        help='A flag to specify if the reads should be mapped  \
                        against human, cat, dog and mouse genomes to find \
                        contaminants. Default = True.')

    parser.add_argument('--paired', '-p', action='store_true', default=False,
                        help='A flag to specify if the fastq file is \
                        paired and interleaved. Default = True.')

    args = parser.parse_args()
    return args

def main():

    args = myparser()

    # Utility functions
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

    # Create the final output file
    if not os.path.isdir(args.output):
        os.mkdir(args.output)

    # Set up logging
    logging.basicConfig(filename=os.path.join(args.output, 'rqc.log'),
                        level=logging.INFO,
                        format='%(asctime)s %(message)s')
    logging.info('Starting quality control workflow.')

    # Assign globals
    abs_fastq = os.path.abspath(args.fastq)

    # Create temporary directory
    rqctempdir = tempfile.mkdtemp()

    # Remove contaminants
    tmp_fc = mk_temp_dir(rqctempdir, 'filter_contaminants')  # make temp dir
    logging.info('Starting contaminant removal')
    data1 = rqcmain.Fastq(path=abs_fastq)
    deconfiles = data1.filter_contaminants(outdir=tmp_fc)
    logging.info(deconfiles)

    # Trim adaptors
    tmp_ta = mk_temp_dir(rqctempdir, 'trim_adaptors')  # make temp. dir.
    logging.info('Starting adaptor trimming')
    data2 = rqcmain.Fastq(os.path.join(tmp_fc, 'clean1.fq.gz'))
    trimfiles = data2.trim_adaptors(outdir=tmp_ta)
    logging.info(trimfiles)

    # Remove vertebrate contaminants
    if args.removevertebrates:
        tmp_rvc = mk_temp_dir(rqctempdir, 'remove_vertebrate_contaminants')
        logging.info('Removing dog, cat, mouse and humanreads')
        # Create new Fastq object
        data3 = rqcmain.Fastq(os.path.join(tmp_ta, 'clean2.fq.gz'))
        if not os.path.isdir(data/dogcatmousehuman/ref):
            pass  # pass until this is on the server with the references
            # rqcmain.build_vertebrate_db(cat=data/cat.fa.gz,
            #                     dog=data/dog.fa.gz,
            #                     human=data/hg19.fa.gz,
            #                     mouse=data/mouse.fa.gz,
            #                     datadir=data/dogcatmousehuman
        rvcfiles = data3.remove_vertebrate_contaminants(outdir=tmp_rvc)
        logging.info(rvcfiles)

    # Merge files
    if args.paired:
        tmp_mr = mk_temp_dir(rqctempdir, 'merge_reads')  # make temp. dir.
        logging.info('Merging read pairs')
        # determine iput file
        if args.removevertebrates:
            infile = os.path.join(tmp_vcr, 'novert.fq.gz')
        else:
            infile = os.path.join(tmp_ta,'clean2.fq.gz')

        data4 = rqcmain.Fastq(os.path.join(tmp_mr, infile))
        mergefiles = data4.merge_reads(outdir=tmp_mr)
        logging.info(mergefiles)

    # For now, copy all files from the tempdir to the output dir
    logging.info("Copying files from temporary directory to ouput directory")
    shutil.copytree(rqctempdir, os.path.join(args.output, "output"))
    shutil.rmtree(rqctempdir)
    logging.info("Completed RQC run")

if __name__ == '__main__':
    main()
