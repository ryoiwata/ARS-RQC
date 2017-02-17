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


def build_vertebrate_db(cat, dog, human, mouse, datadir):
    """Builds a bbmap.sh database for mapping reads to masked versions of
    the cat, dog, human and mouse genome. Returns the path of the database"""
    try:
        parameters = ['bbmap.sh', 'build=1', 'k=14', 'usemodulo',
                      'ref=' +
                      os.path.abspath(cat) + ',' +
                      os.path.abspath(dog) + ',' +
                      os.path.abspath(human) + ',' +
                      os.path.abspath(mouse),
                      'path=' + os.path.join(datadir, 'dogcatmousehuman')]
        subprocess.run(parameters)
        return os.path.join(datadir, 'dogcatmousehuman')
    except RuntimeError:
        print("Couldn't build database of vertebrate contaminants using bbmap")


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
        """Calls bbduk to perform adapter removal"""
        try:
            bbtoolsdict = self.parse_params()
            parameters = ['bbduk.sh',
                          'in=' + self.abspath,
                          'out=' + os.path.join(outdir, 'clean1.fq.gz'),
                          'outduk=' + os.path.join(outdir, 'kmerStats1.txt'),
                          'stats=' + os.path.join(outdir, 'scaffoldStats1.txt')]
            parameters.extend(bbtoolsdict['filter_contaminants'])
            p1 = subprocess.run(parameters, check=True, stderr=subprocess.PIPE)
            self.metadata['filter_contaminants'] = os.walk(outdir)
        except RuntimeError:
            print("could not perform contaminant filtering with bbduk")

    def trim_adaptors(self, outdir):
        """Calls bbduk to remove contaminant sequences"""
        try:
            bbtoolsdict = self.parse_params()
            parameters = ['bbduk.sh', 'in=' + self.abspath,
                          'out=' + os.path.join(outdir, 'clean2.fq.gz'),
                          'stats=' + os.path.join(outdir, 'scaffoldStats2.txt')]
            parameters.extend(bbtoolsdict['trim_adaptors'])
            p2 = subprocess.run(parameters)
            self.metadata['trim_adaptors'] = os.walk(outdir)
        except RuntimeError:
            print("could not perform adaptor removal with bbduk")

    def merge_reads(self, outdir):
        """merge reads to generate insert size histogram and \
        optionally save"""
        try:
            bbtoolsdict = self.parse_params()
            parameters = ['bbmerge.sh', 'in1=' + self.abspath,
                          'ihist=' + os.path.join(outdir, 'merge_histogram.txt'),
                          'outc=' + os.path.join(outdir, 'cardinality.txt'),
                          'outm=' + os.path.join(outdir, 'merged.fq.gz'),
                          'outu=' + os.path.join(outdir, 'unmerged.fq.gz')]
            parameters.extend(bbtoolsdict['merge_reads'])
            subprocess.run(parameters)
            self.metadata['merge_reads'] = os.walk(outdir)
        except RuntimeError:
            print("could not perform read merging with bbmerge")

    def remove_vertebrate_contaminants(self, outdir):
        """maps reads to repeat-masked human, dog, cat and mouse genomes
        to remove contaminants"""
        try:
            bbtoolsdict = self.parse_params()
            parameters = ['bbmap.sh',
                          'in=' + self.abspath,
                          'outu=' + os.path.join(outdir, 'novert.fq.gz')]
            parameters.extend(bbtoolsdict['remove_vertebrate_contaminants'])
            suborocess.run(parameters)
            self.metadata['remove_vertebrate_contaminants'] = os.walk(outdir)
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
            subprocess.run(['sortbyname.sh',
                            "in=" + self.abspath,
                            "out=" + tfile])
            shutil.move(tfile, os.path.join(root, uc))
        except RuntimeError:
            print("could not reorder fastq file by name")
        finally:
            shutil.rmtree(temp_ordered_dir)

def main():
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
                        help='A flag to specify if the reads should be mapped  \
                        against human, cat, dog and mouse genomes to find \
                        contaminants. Default = True.')

    parser.add_argument('--paired', '-p', action='store_true',
                        help='A flag to specify if the fastq file is \
                        paired and interleaved. Default = True.')

    parser.add_argument('--keepmerged', '-km', action='store_true',
                        help='Keep fastq files with merged and unmerged reads. \
                        Defualt = True. and unmerged reads will be output')

    args = parser.parse_args()

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

    # Validate options
    if args.fastqpair and args.fastq is not None:
        raise Exception("Please pass only a paired or unpaired fastq or \
                        fastq.gz file, not both.")
        logging.error('Please pass only a paired or unpaired fastq or fastq.gz \
                      file, not both.')

    # Assign globals
    abs_fastq = os.path.abspath(args.fastq)
    abs_datadir = os.path.abspath('data')

    # Create temporary directory
    rqctempdir = tempfile.TemporaryDirectory()

    # Remove contaminants
    tmp_cf = mk_temp_dir(rqctempdir, 'contaminant_filtering')  # make temp dir
    logging.info('Starting contaminant removal')
    data1 = Fastq(path=abs_fastq)
    deconfiles = data1.filter_contaminants(tmp_cf)

    # Trim adaptors
    tmp_ta = mk_temp_dir(rqctempdir, 'trim_adaptors')  # make temp. dir.
    logging.info('Starting adaptor trimming')
    data2 = Fastq(os.path.join(tmp_cf, 'clean1.fq.gz'))
    trimfiles = data.trim_adaptors(path=tmp_ta)

    if args.removevertebrates:
        tmp_rv = mk_temp_dir(rqctempdir, 'vertebrate_contaminant_removal')
        logging.info('Removing dog, cat, mouse and humanreads')
        # Create new Fastq object
        data3 = Fastq(os.path.join(tmp_rv, 'clean2.fq.gz'))
        if not os.path.isdir(data/dogcatmousehuman/ref):
            pass  # pass until this is on the server with the references
            build_vertebrate_db(cat=data/cat.fa.gz,
                                dog=data/dog.fa.gz,
                                human=data/hg19.fa.gz,
                                mouse=data/mouse.fa.gz,
                                datadir=data/dogcatmousehuman)

        trimfiles = data.remove_vertebrate_contaminants(path=tmp_rv)

    if args.paired:
        # merge files
        tmp_mr = mk_temp_dir(rqctempdir, 'merge_reads')  # make temp. dir.
        logging.info('Merging reads pairs')

        # Create new Fastq object
        if args.removevertebrates:
            infile = 'novert.fq.gz'
        else:
            infile = 'clean2.fq.gz'
        data4 = Fastq(os.path.join(tmp_mr, infile))
        mergefiles = data.merge_reads(path=infile)

    # For now copy all files from the tempdir to the output dir
    shutil.copytree(rqctmpdir, args.output)


if __name__ == '__main__':
    main()
