#!/bin/bash
#SBATCH -p short #name of the partition (queue) you are submitting to
#SBATCH -N 1 #number of nodes in this job
#SBATCH -n 40 #number of cores/tasks in this job,
#SBATCH -t 02:00:00 #time allocated for this job hours:mins:seconds
#SBATCH --mail-user=emailAddress #enter your email address to receive emails
#SBATCH --array=1-$1
# usage ceresjobsub.sh  filenum <indir> <outdir>


file=`ls "$2" | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`
fname=`basename $file`
rqcparser.py --fastq "$file" --output "$3"/"$fname" -m -p
