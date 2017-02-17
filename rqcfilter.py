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


# def build_vertebrate_db(cat, dog, human, mouse, datadir):
#     """Builds a bbmap.sh database for mapping reads to masked versions of
#     the cat, dog, human and mouse genome. Returns the path of the database"""
#     try:
#         parameters = ['bbmap.sh', 'build=1', 'k=14', 'usemodulo',
#                       'ref=' +
#                       os.path.abspath(cat) + ',' +
#                       os.path.abspath(dog) + ',' +
#                       os.path.abspath(human) + ',' +
#                       os.path.abspath(mouse),
#                       'path=' + os.path.join(datadir, 'dogcatmousehuman')]
#         subprocess.run(parameters)
#         return os.path.join(datadir, 'dogcatmousehuman')
#     except RuntimeError:
#         print("Couldn't build database of vertebrate contaminants using bbmap")
#
#
# class Fastq():
#
#     @classmethod
#     def parse_params(self):
#         """converts parameters from a json file into a dictionary formatted
#         for use by bbtools"""
#         try:
#             with open("data/parameters.json", 'r') as p:
#                 fastq_parameters = json.load(p)
#             bbtoolsdict = {}
#             for key in fastq_parameters["rqcfilter"]:
#                 bbtoolsdict[key] = []
#                 for key2 in fastq_parameters["rqcfilter"][key]:
#                     if isinstance(fastq_parameters["rqcfilter"][key][key2],
#                                   str):
#                         bbtoolsdict[key].append(str(key2) + "=" +
#                                                 str(fastq_parameters
#                                                 ["rqcfilter"][key][key2]))
#                     elif isinstance(fastq_parameters["rqcfilter"][key][key2],
#                                     list):
#                         bbtoolsdict[key].append(str(key2) + "=" +
#                                                 ",".join(fastq_parameters
#                                                          ["rqcfilter"]
#                                                          [key][key2]))
#             return bbtoolsdict
#         except RuntimeError:
#                 print("Could not load and parse the parameters.json file \
#                       correctly")
#
#     def __init__(self, path):
#         self.abspath = os.path.abspath(path)
#         self.filepath, self.filename = os.path.split(os.path.abspath(path))
#         self.metadata = {}
#
#     def __repr__(self):
#         return 'Fastq Class object :' + self.filename
#
#     def filter_contaminants(self, outdir):
#         """Calls bbduk to perform adapter removal"""
#         try:
#             bbtoolsdict = self.parse_params()
#             parameters = ['bbduk.sh',
#                           'in=' + self.abspath,
#                           'out=' + os.path.join(outdir, 'clean1.fq.gz'),
#                           'outduk=' + os.path.join(outdir, 'kmerStats1.txt'),
#                           'stats=' + os.path.join(outdir, 'scaffoldStats1.txt')]
#             parameters.extend(bbtoolsdict['filter_contaminants'])
#             p1 = subprocess.run(parameters, check=True, stderr=subprocess.PIPE)
#             self.metadata['filter_contaminants'] = list(os.walk(outdir))
#             return p1.stderr.decode('utf-8')
#         except RuntimeError:
#             print("could not perform contaminant filtering with bbduk")
#
#     def trim_adaptors(self, outdir):
#         """Calls bbduk to remove contaminant sequences"""
#         try:
#             bbtoolsdict = self.parse_params()
#             parameters = ['bbduk.sh', 'in=' + self.abspath,
#                           'out=' + os.path.join(outdir, 'clean2.fq.gz'),
#                           'stats=' + os.path.join(outdir, 'scaffoldStats2.txt')]
#             parameters.extend(bbtoolsdict['trim_adaptors'])
#             p2 = subprocess.run(parameters, check=True, stderr=subprocess.PIPE)
#             self.metadata['trim_adaptors'] = list(os.walk(outdir))
#             return p2.stderr.decode('utf-8')
#         except RuntimeError:
#             print("could not perform adaptor removal with bbduk")
#
#     def merge_reads(self, outdir):
#         """merge reads to generate insert size histogram and \
#         optionally save"""
#         try:
#             bbtoolsdict = self.parse_params()
#             parameters = ['bbmerge.sh', 'in1=' + self.abspath,
#                           'ihist=' + os.path.join(outdir, 'merge_histogram.txt'),
#                           'outc=' + os.path.join(outdir, 'cardinality.txt'),
#                           'outm=' + os.path.join(outdir, 'merged.fq.gz'),
#                           'outu=' + os.path.join(outdir, 'unmerged.fq.gz')]
#             parameters.extend(bbtoolsdict['merge_reads'])
#             p3 = subprocess.run(parameters, check=True, stderr=subprocess.PIPE)
#             self.metadata['merge_reads'] = list(os.walk(outdir))
#             return p3.stderr.decode('utf-8')
#         except RuntimeError:
#             print("could not perform read merging with bbmerge")
#
#     def remove_vertebrate_contaminants(self, outdir):
#         """maps reads to repeat-masked human, dog, cat and mouse genomes
#         to remove contaminants"""
#         try:
#             bbtoolsdict = self.parse_params()
#             parameters = ['bbmap.sh',
#                           'in=' + self.abspath,
#                           'outu=' + os.path.join(outdir, 'novert.fq.gz')]
#             parameters.extend(bbtoolsdict['remove_vertebrate_contaminants'])
#             p4 = suborocess.run(parameters, check=True, stderr=subprocess.PIPE)
#             self.metadata['remove_vertebrate_contaminants'] = list(os.walk(outdir))
#             return p4.stderr.decode('utf-8')
#         except RuntimeError:
#             print("Could not perform vertebrate conaminant removal with bbmap")
#
#     def sortbyname(self):
#         """Sorts a fastq file by read names, outputs uncompressed fastq"""
#         try:
#             temp_ordered_dir = tempfile.mkdtemp()
#             root, base = os.path.split(self.abspath)
#             if base.endswith(".gz"):
#                 uc = base.rsplit(".", 1)[0]
#             else:
#                 uc = base
#             tfile = os.path.join(temp_ordered_dir, 'sorted.fq')
#             subprocess.run(['sortbyname.sh',
#                             "in=" + self.abspath,
#                             "out=" + tfile])
#             shutil.move(tfile, os.path.join(root, uc))
#         except RuntimeError:
#             print("could not reorder fastq file by name")
#         finally:
#             shutil.rmtree(temp_ordered_dir)


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
