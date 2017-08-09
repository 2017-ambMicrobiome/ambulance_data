#!/bin/bash
#$ -cwd
#$ -j n
#$ -l h_rt=96:00:00
#$ -m bea
#$ -M eba2001@med.cornell.edu
#$ -N ambulance_metaphlan
#$ -pe smp 5
#$ -l os=rhel6.3
#$ -l h_vmem=1G

rsync -a /zenodotus/masonlab/pathomap_scratch/ebrahim/Ambulance_Study/clean_up/*fastq.gz $TMPDIR
mkdir $TMPDIR/metaphlan_out
cd $TMPDIR
for file in $(find $TMPDIR/metaphlan_out *fastq.gz)
do
zcat $file | /home/ebrahim/anaconda/bin/python /home/ebrahim/bin/metaphlan2/metaphlan2.py --bowtie2db /home/ebrahim/bin/metaphlan2/db_v20/mpa_v20_m200 --bowtie2_exe /home/ebrahim/bin/bowtie2-2.1.0/bowtie2 --input_type fastq --mpa_pkl /home/ebrahim/bin/metaphlan2/db_v20/mpa_v20_m200.pkl --nproc $NSLOTS --bowtie2out $file.bt2.out > $TMPDIR/metaphlan_out/$file.out

done

rsync -a $TMPDIR/metaphlan0_out /zenodotus/masonlab/pathomap_scratch/ebrahim/Ambulance_Study/