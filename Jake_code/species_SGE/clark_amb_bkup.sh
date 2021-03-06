#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_rt=96:00:00
#$ -m bea
#$ -M harryjreed1985@gmail.com
#$ -N clark_d
#$ -pe smp 16
#$ -l os=rhel6.3
#$ -l vf=10G
#$ -l h_vmem=16G
#$ -l zenodotus=true

# Path variables
pathFromZeno=/zenodotus/masonlab/rotcollab_scratch/nbo5/trim_amb_samples/trim_amb_adapt
pathToZeno=/zenodotus/masonlab/rotcollab_scratch/har2011/ambulance
pathToCLARK=/home/har2011/bin/CLARKSCV1.2.2-b
pathToDB=/home/har2011/bin/clark_db

# Sync data to tmpdir
cd $TMPDIR

echo -e "\nclark_ambulance_debug.sh started on ambulance samples at `date`.\nSyncing from: ${pathZeno} \nTo: ${TMPDIR}/clark01_work\n"

mkdir clark01_out
mkdir clark01_work
mkdir clark01_abundance_out
mkdir clark

cd ${pathFromZeno}
find . -maxdepth 1 *.gz -exec ls -lh {} \+ | egrep "./AE|./AS|./AW" > ref
rsync -av $( cat ref | grep -o './.*' | head -n 4 )  ${TMPDIR}/clark01_work
rm ref
rsync -av ${pathToCLARK}/* ${TMPDIR}/clark/

# Unzipping files
echo -e "\nStarted to gunzip files at `date`.\n"
cd $TMPDIR/clark01_work
gunzip *
echo -e "Finished gunzip files at `date`.\n"

# Setting targets from database
cd $TMPDIR/clark
	./set_targets.sh ${pathToDB} bacteria --species

# Classifying metagenome PE reads
echo -e "\nMetagenome classification began at `date`.\n"
cd $TMPDIR/clark
for file in $( find ${TMPDIR}/clark01_work/* | grep -o 'A....._' | awk 'NR%2==0' )
do 
	./classify_metagenome.sh -P ${TMPDIR}/clark01_work/${file}pair01.fastq.out.clean ${TMPDIR}/clark01_work/${file}pair02.fastq.out.clean -R ${TMPDIR}/clark01_out/${file}clarkOut -m 1 -n 16 --ldm
done
echo -e "\nMetagenome classification completed at `date`.\n"

# Abundance estimation of output files with high-confidence
echo -e "Abundance estimation began at `date`.\n"
cd $TMPDIR/clark
for clarkFile in $( find ${TMPDIR}/clark01_out/*.csv )
do
	./estimate_abundance.sh -F ${clarkFile} -D ${pathToDB} -c 0.90 -g 0.03 -a 0.001 > ${clarkFile}_abundance
done

# Mv files and rename to clark_01_abundance_out
mv ${TMPDIR}/clark01_out/*_abundance ${TMPDIR}/clark01_abundance_out/
cd ${TMPDIR}/clark01_abundance_out/
for rename in $( ls | egrep "AE|AS|AW" | grep -o 'A.....' ) 
do 
	mv ${rename}_clarkOut.csv_abundance ./${rename}_abundance.csv 
	echo -e "\nRenaming file ${rename}."
done
echo -e "\nAbundance estimation and renaming finished at `date`.\n"

# Zipping clark01_out files
echo -e "Zipping files started at `date`.\n"
cd $TMPDIR/clark01_out
gzip *
echo -e "\nFinished zipping files at `date`.\n"

# Sync data back to $myamdata
echo -e "Syncing files back to zeno.\n"
rsync -av ${TMPDIR}/clark01_out/* ${pathToZeno}/CLARK_out
rsync -av ${TMPDIR}/clark01_abundance_out/* ${pathToZeno}/CLARK_abundance_out

echo -e "\nFinished running clark_ambulance_debug.sh at `date`."
