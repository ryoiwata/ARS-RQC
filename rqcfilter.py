#!/usr/env/python3
# rqcfilter.py - A sequence quality control and metadata collection workflow

import argparse
import os
import subprocess
import logging
import tempfile
import json

def create_parameters(dict):
    """creates a dictionary containing formatted parameters in list format"""


def is_gzipped(path):
    """determine if fastq file is gzipped by its extension"""
    extension=os.path.splitext(path)[1]
    if extension == 'gz':
        return True
    elif extension == 'fastq' or 'fq':
        return False
    else:
        raise Exception("Only files ending in .fastq.gz, fq.gz, .fastq \
        or .fq can be used")

def modulo(f):
    name = f.readline()
    seq = f.readline()
    assert type(seq) is StringType, 'sequence is not a string: {}'.format(seq)
    return len(seq) % 5 == 0

def is_fastq_modulo_zero(fastq):
    """check the length of a compressed or uncompressed fastq file"""
    if is_gzipped(fastq):
        import gzip
        with gzip.open(fastq, 'rb') as f:
            return modulo(f)
    else:
        with open(fastq, 'r') as f:
            return modulo(f)

def trim_adaptors(fastq, tmpdir, config):
    """Checks sequence length and runs bbduk for adapter removal"""

    parameters = ['bbduk.sh','in=' + fastq,
                'out='  + tmpdir + '/clean.fq']
                "ref=data/Illumina.artifacts.2013.12.no_DNA_RNA_spikeins.fa, \
                data/DNA_spikeins.artifacts.2012.10.fa, \
                data/phix174_ill.ref.fa,data/lambda.fa.gz, \
                data/pJET1.2.fasta"]

    parameters.append(trim_paramters)
    if is_fastq_modulo_zero(fastq) == False:
        parameters = parameters.append('ftm=5')
    p1 = subprocess.run(parameters, stderr=subprocess.PIPE)
    s1 = subprocess.communitate(p1)
    return s1

def contaminant_filtering(fastq, tmpdir):
    """removes contaminant sequences removal"""
    parameters = ['bbduk.sh', 'in=' + tmpdir + '/clean.fq',
    'out='  + tmpdir + '/clean_decon.fq']
    parameters.append(decon_paramters, stderr=subprocess.PIPE)
    p2 = subprocess.run(parameters)
    s2 = subprocess.communitate(p2)
    return s2
def mergereads(fastq,tmpdir):
    """merge reads to generate insert size histogram and optionally save """
    parameters = ['bbmerge.sh', 'in1=' + tmpdir + '/clean_decon.fq', \
    'ihist=' + tmpdir + '/ihist_merge.txt', 'outc=' + tmpdir+ '/cardinality.txt']
    parameters = parameters.append()
    parameters.append(bbmerge_parameters)
    if args.keepmerged == True:
        parameters.append(['outm=' +tmpdir + '/clean_decon_m.fq.gz]', \
        'outu=' + tmpdir + '/clean_decon_u.fq.gz'])
    subprocess.run(parameters)

def build_vertebrate_db(cat,dog,human,mouse):
    """builds a bbmap.sh database for mapping reads to masked versions of \
    the cat, dog, human and mouse genome"""
    parameters = ['bbmap.sh', 'build=1', 'k=14', 'usemodulo', 'ref=' + os.path.abspath(cat) + \
    ',' + os.path.abspath(cat) + ',' + 'os.path.abspath(human)' + ',' + \
    os.path.abspath(mouse), 'path=' + os.path.join(abs_datadir,'dogcatmousehuman')]
    subprocess.run(parameters)

def vertebrate_contaminant_removal(fastq,tmpdir):
    """maps reads to repeat-masked human, dog, cat and mouse genomes \
    to remove contaminants"""
    parameters = ['bbmap.sh', 'build=1','in=' + fastq, 'outu=' + tmpdir + \
    'clean_decon_novert.fq.gz']
    parameters.append(bbmap_parameters)
    supbrosess.run(parameters)


def main:
    parser = argparse.ArgumentParser(description='rqcfilter.py - A sequence quality \
                        control and metadata collection workflow.')

    parser.add_argument('--fastq', '-f', type=str, required=True,
                        help='A .fastq, .fq, .fastq.gz or .fq.gz file.')

    parser.add_argument('--parameters.json', '-m', type=argparse.FileType('r'),
                        help='A yaml formated metadata file. See documentation.')

    parser.add_argument('--qcdata', '-qc', type=str, default='qcdata.json',
                        help='The name of the file containing QC data in json format')

    parser.add_argument('--output', '-o', type=str, default='rqcout',
                        help='the output directory')

    parser.add_argument('--datadir', '-d', type=str, default='rqcout',
                        help='the output directory')

    parser.add_argument('--paired', '-p', action='store_true',
                        help='A flag to specify if the fastq file is \
                        paired and interleaved. Default = True.')

    parser.add_argument('--keepmerged', '-km', action='store_true',
                        help='Keep fastq files with merged and unmerged reads. Defualt = True. \
                        and unmerged reads will be output')

    parser.add_argument('--threads', '-t', type=int, default=16,
                        help='The Number of threads for multithreaded components to use.')

    args = parser.parse_args()

    if not os.path.isdir(args.output):
        os.mkdir(args.output)

    # set up logging
    logging.basicConfig(filename=str(args.output) + '/rqc.log', level=logging.INFO, \
                        format='%(asctime)s %(message)s')
    logging.info('Starting quality control workflow.')

    # Validate options
    if args.fastqpair and args.fastq != None:
        raise Exception("Please pass only a paired or unpaired fastq or fastq.gz \
            file, not both.")
        logging.error('Please pass only a paired or unpaired fastq or fastq.gz file,\
            not both.')

    # Assign globals
    abs_fastq = os.path.abspath(args.fastq)
    abs_datadir = os.path.abspath('data')
    os.mkdir(os.path.abspath(args.output))

if __name__ == '__main__':
    main()
