
#!/bin/bash
## config.sh, a script to set up the qc workflow on CERES
# Adam R. Rivers, USDA - ARS - GBRU
# 01/24/2017

# # Set up bbtools - no longer needed now that module has been installed
# mkdir RQC
# cd RQC
# mkdir bin
# # load java sdk
# module load java_sdk/64/1.7.0_79
# export JAVA_HOME="/software/apps/java_sdk/64/1.7.0_79/bin/java"
# wget "https://sourceforge.net/projects/bbmap/files/latest/download" -0
# tar -xvf download
# rm -r download
# mv bbmap bin/
# cd bin/bbmap/jni
# # create c++ binaries
# make -f makefile.linux
# cd ../../

module load bbmap

# Download pigz for parallel (de)compression
wget "http://zlib.net/pigz/pigz-2.3.4.tar.gz"
tar -xvf pigz-2.3.4.tar.gz
cd pigz
make
cd ../


bbduk.sh maq=5 trimq=10 qtrim=f overwrite=true maxns=0 minlen=25 minlenfraction=0.333 k=31 hdist=1 pigz=t unpigz=t zl=6 cf=t barcodefilter=crash ow=true rqc=hashmap outduk=test/kmerStats1.txt stats=test/scaffoldStats1.txt loglog ref=data/Illumina.artifacts.2013.12.no_DNA_RNA_spikeins.fa,data/DNA_spikeins.artifacts.2012.10.fa,data/phix174_ill.ref.fa,data/lambda.fa.gz,data/pJET1.2.fasta in=data/sample.fq.gz 2> sample_trim_adaptors.txt
