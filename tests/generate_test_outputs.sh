#!/bin/bash
bbduk.sh \
maq=5 \
trimq=10 \
qtrim=f \
overwrite=true \
maxns=0 \
minlen=25 \
minlenfraction=0.333 \
k=31 \
hdist=1 \
pigz=t \
unpigz=t \zl=6 \
cf=t \
barcodefilter=crash \
ow=true \
rqc=hashmap \
outduk=test/kmerStats1.txt \
stats=test/scaffoldStats1.txt \
loglog \
ref=../data/Illumina.artifacts.2013.12.no_DNA_RNA_spikeins.fa,../data/DNA_spikeins.artifacts.2012.10.fa,../data/phix174_ill.ref.fa,../data/lambda.fa.gz,../data/pJET1.2.fasta in=../data/reads100.fq.gz out=sample.trimmed.fq.gz \
threads=1 \
interleaved=t \
-Xmx4g 2> sample_trim_adaptors.txt
