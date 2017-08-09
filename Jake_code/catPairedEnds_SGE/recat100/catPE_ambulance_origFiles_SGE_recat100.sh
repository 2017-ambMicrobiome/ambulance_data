#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_rt=96:00:00
#$ -m bea
#$ -M harryjreed1985@gmail.com
#$ -N recat100
#$ -pe smp 8
#$ -l os=rhel6.3
#$ -l vf=1G
#$ -l h_vmem=2G
#$ -l zenodotus=true

# Setting Variables
pathToSync=/zenodotus/masonlab/pathomap_scratch/ebrahim/Ambulance_Study/hudson_alpha/FASTQ
pathZenoSync=/zenodotus/masonlab/rotcollab_scratch/har2011/ambulance/FASTQ_hudson_alpha/pairedEnds_fromOriginal/one_1

cd $TMPDIR
mkdir finished

# Syncing Data
echo -e "Started syncing data at `date`\n"
rsync -av $( find $pathToSync/* | egrep "C6R85ANXX|HADEY" | grep -v "BC_X" )  $TMPDIR

# Syncing refFiles
rsync -av /home/har2011/SGE_scripts/catPairedEnds/recat100/RefFile100 $TMPDIR

# catPair01_ambulance
echo -e "Started to recat100 at `date`\n"
bash RefFile100
mv A* finished/
echo -e "Finished recatting100 at `date`\n" 

# syncing back to zenodotus
rsync -av $TMPDIR/finished/* $pathZenoSync 
