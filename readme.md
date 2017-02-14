This repository contains a quality control/ quality assurance system for sequencing data.

Inputs
1. demultiplexed fastq files
2. metadata
3. runs BBtools RQC pipeline
4. collects metadata about the quality of the reads
5. runs analytics and stores the data
  * kmer spectrum
  * GC divergence
  * GATC by position
  * insert size
  * sequence GC distribution
