#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_rt=96:00:00
#$ -m bea
#$ -M harryjreed1985@gmail.com
#$ -N clark_d
#$ -pe smp 4
#$ -l os=rhel6.3
#$ -l vf=35G
#$ -l h_vmem=55G
#$ -l zenodotus=true

# Path variables
pathFromZeno=/zenodotus/masonlab/rotcollab_scratch/har2011/ambulance/FASTQ_hudson_alpha/repair_PEs/pairs
pathToZeno=/zenodotus/masonlab/rotcollab_scratch/har2011/ambulance
pathToSoft=/home/har2011/bin
pathToDB=/home/har2011/bin/clark_db

# Sync data to tmpdir
cd $TMPDIR

echo -e "\nclark_ambulance_debug.sh started on ambulance samples at `date`.\nSyncing from: ${pathZeno} \nTo: ${TMPDIR}/clark01_work\n"

mkdir clark01_out
mkdir clark01_work
mkdir clark01_abundance_out
mkdir soft

cd ${pathFromZeno}
rsync -av $( ls | egrep "AE|AS|AW" | grep -o "A.*" | head -n 4 )  ${TMPDIR}/clark01_work
rsync -av ${pathToSoft}/CLARKSCV1.2.2-b ${TMPDIR}/soft
rsync -av ${pathToSoft}/parallel ${TMPDIR}/soft

# Unzipping files
echo -e "\nStarted to gunzip files at `date`.\n"
cd $TMPDIR/clark01_work
ls | ${TMPDIR}/soft/parallel --gnu -j $NSLOTS gunzip
echo -e "Finished gunzip files at `date`.\n"

# Setting targets from database
cd $TMPDIR/soft/CLARKSCV1.2.2-b 
	./set_targets.sh ${pathToDB} bacteria --species

# Classifying metagenome PE reads
echo -e "\nMetagenome classification began at `date`.\n"
cd $TMPDIR/soft/CLARKSCV1.2.2-b 
for file in $( find ${TMPDIR}/clark01_work/* | grep -o 'A....._' | awk 'NR%2==0' )
do 
	./classify_metagenome.sh -P ${TMPDIR}/clark01_work/${file}pairs_R1.fastq ${TMPDIR}/clark01_work/${file}pairs_R2.fastq -R ${TMPDIR}/clark01_out/${file}clarkOut -m 1 -n $NSLOTS --ldm
done
echo -e "\nMetagenome classification completed at `date`.\n"

# Abundance estimation of output files with high-confidence
echo -e "Abundance estimation began at `date`.\n"
cd $TMPDIR/soft
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
find . -name '*.csv' | sort > gzRef 
${TMPDIR}/soft/parallel --gnu -j $NSLOTS gzip :::: gzRef
rm gzRef
echo -e "\nFinished zipping files at `date`.\n"

# Sync data back to $myamdata
echo -e "Syncing files back to zeno.\n"
rsync -av ${TMPDIR}/clark01_out/* ${pathToZeno}/CLARK_out
rsync -av ${TMPDIR}/clark01_abundance_out/* ${pathToZeno}/CLARK_abundance_out

echo -e "\nFinished running clark_ambulance_debug.sh at `date`."
