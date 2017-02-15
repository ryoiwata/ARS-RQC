#!/usr/env/python3
# rqcfilter.py - A sequence quality control and metadata collection workflow

import argparse
import os
import subprocess
import logging
import tempfile
import json


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

    def trim_adaptors(self, outdir):
        """Calls bbduk to perform adapter removal"""
        try:
            bbtoolsdict = self.parse_params()
            parameters = ['bbduk.sh', 'in=' + self.abspath, 'out=' +
                          os.path.join(outdir, 'clean.fq.gz'),
                          "outduk=" + os.path.join(outdir, )]
            parameters.extend(bbtoolsdict["trim_adaptors"])
            p1 = subprocess.run(parameters, check=True, stderr=subprocess.PIPE)
            # print(p1.stderr)
            # return p1
        except RuntimeError:
            print("could not perform adaptor removal")

    def contaminant_filtering(fastq, tmpdir):
        """Calls bbduk to remove contaminant sequences"""
        parameters = ['bbduk.sh', 'in=' + tmpdir + '/clean.fq',
                      'out=' + tmpdir + '/clean_decon.fq']
        parameters.append(decon_paramters, stderr=subprocess.PIPE)
        p2 = subprocess.run(parameters)
        s2 = subprocess.communitate(p2)
        return s2

    def mergereads(fastq, tmpdir):
        """merge reads to generate insert size histogram and
        optionally save"""
        parameters = ['bbmerge.sh', 'in1=' + tmpdir + '/clean_decon.fq',
                      'ihist=' + tmpdir + '/ihist_merge.txt', 'outc=' +
                      tmpdir + '/cardinality.txt']
        parameters = parameters.append()
        parameters.append(bbmerge_parameters)
        if args.keepmerged is True:
            parameters.append(['outm=' + tmpdir + '/clean_decon_m.fq.gz]',
                               'outu=' + tmpdir + '/clean_decon_u.fq.gz'])
        subprocess.run(parameters)

    def build_vertebrate_db(cat, dog, human, mouse):
        """builds a bbmap.sh database for mapping reads to masked versions of
        the cat, dog, human and mouse genome"""
        parameters = ['bbmap.sh', 'build=1', 'k=14', 'usemodulo',
                      'ref=' + os.path.abspath(cat) + ',' +
                      os.path.abspath(cat) +
                      ',' + 'os.path.abspath(human)' + ',' +
                      os.path.abspath(mouse),
                      'path=' + os.path.join(abs_datadir, 'dogcatmousehuman')]
        subprocess.run(parameters)

    def vertebrate_contaminant_removal(fastq, tmpdir):
        """maps reads to repeat-masked human, dog, cat and mouse genomes
        to remove contaminants"""
        parameters = ['bbmap.sh', 'build=1', 'in=' + fastq, 'outu=' + tmpdir +
                      'clean_decon_novert.fq.gz']
        parameters.append(bbmap_parameters)
        supbrosess.run(parameters)


def main():
    parser = argparse.ArgumentParser(description='rqcfilter.py - \
                                     A sequence quality control and metadata \
                                     collection workflow.')

    parser.add_argument('--fastq', '-f', type=str, required=True,
                        help='A .fastq, .fq, .fastq.gz or .fq.gz file.')

    parser.add_argument('--parameters.json', '-m', type=argparse.FileType('r'),
                        help='A yaml formated metadata file. \
                        See documentation.')

    parser.add_argument('--qcdata', '-qc', type=str, default='qcdata.json',
                        help='The name of the file containing QC data in \
                        json format')

    parser.add_argument('--output', '-o', type=str, default='rqcout',
                        help='the output directory')

    parser.add_argument('--datadir', '-d', type=str, default='rqcout',
                        help='the output directory')

    parser.add_argument('--paired', '-p', action='store_true',
                        help='A flag to specify if the fastq file is \
                        paired and interleaved. Default = True.')

    parser.add_argument('--keepmerged', '-km', action='store_true',
                        help='Keep fastq files with merged and unmerged reads. \
                        Defualt = True. and unmerged reads will be output')

    parser.add_argument('--threads', '-t', type=int, default=16,
                        help='The Number of threads for multithreaded \
                        components to use.')

    args = parser.parse_args()

    if not os.path.isdir(args.output):
        os.mkdir(args.output)

    # set up logging
    logging.basicConfig(filename=str(args.output) + '/rqc.log',
                        level=logging.INFO,
                        format='%(asctime)s %(message)s')
    logging.info('Starting quality control workflow.')

    # Validate options
    if args.fastqpair and args.fastq is not None:
        raise Exception("Please pass only a paired or unpaired fastq or \
                        fastq.gz file, not both.")
        logging.error('Please pass only a paired or unpaired fastq or fastq.gz \
                      file, not both.')

    # Assign globals
    abs_fastq = os.path.abspath(args.fastq)
    abs_datadir = os.path.abspath('data')
    os.mkdir(os.path.abspath(args.output))


if __name__ == '__main__':
    main()
