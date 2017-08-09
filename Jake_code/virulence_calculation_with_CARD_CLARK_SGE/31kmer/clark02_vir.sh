#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_rt=96:00:00
#$ -m bea
#$ -M harryjreed1985@gmail.com
#$ -N clark02_vir
#$ -pe smp 16
#$ -l os=rhel6.3
#$ -l vf=10G
#$ -l h_vmem=15G
#$ -l zenodotus=true

# Path variables
pathFromZeno=/zenodotus/masonlab/rotcollab_scratch/har2011/ambulance/FASTQ_hudson_alpha/repair_PEs/pairs
pathToZeno=/zenodotus/masonlab/rotcollab_scratch/har2011/ambulance
pathToSoft=/home/har2011/software
pathToDB=/zenodotus/masonlab/rotcollab_scratch/har2011/DB_files/clark_virulence_db

# Sync data to tmpdir
cd $TMPDIR

echo -e "\nclark01.sh started on ambulance samples at `date`.\nSyncing from: ${pathZeno} \nTo: ${TMPDIR}/clark01_work\n"

mkdir out
mkdir work
mkdir soft
mkdir soft/clark_virulence_db

# 200-400
rsync -av $( find ${pathFromZeno}/*.gz | sort | head -n 400 | tail -n 200 )  ${TMPDIR}/work
rsync -a ${pathToSoft}/CLARKSCV1.2.2-b ${TMPDIR}/soft
rsync -av ${pathToSoft}/parallel ${TMPDIR}/soft

# Sync only necessary DB files
rsync -a ${pathToDB}/* ${TMPDIR}/soft/clark_virulence_db

# Unzipping files
echo -e "\nStarted to gunzip files at `date`.\n"
cd $TMPDIR/work
ls | ${TMPDIR}/soft/parallel --gnu -j 16 gunzip
echo -e "Finished gunzip files at `date`.\n"

# Classifying metagenome PE reads
cd $TMPDIR/soft/CLARKSCV1.2.2-b/exe
for file in $( find ${TMPDIR}/work/* | grep -o 'A....._' | awk 'NR%2==0' )
do
        ./CLARK -T ${TMPDIR}/soft/clark_virulence_db/targets_finished.txt -D ${TMPDIR}/soft/clark_virulence_db/db_files -P ${TMPDIR}/work/${file}pairs_R1.fastq ${TMPDIR}/work/${file}pairs_R2.fastq -R ${TMPDIR}/out/${file}clarkOut -m 0 -n 16 --ldm -k 31
done

# Zipping clark01_out files
echo -e "Zipping files started at `date`.\n"
cd $TMPDIR/out
find ./*.csv | sort > gzRef 
${TMPDIR}/soft/parallel --gnu -j 16 gzip :::: gzRef
rm gzRef
echo -e "\nFinished zipping files at `date`.\n"

# Sync data back to $myamdata
echo -e "Syncing files back to zeno.\n"
rsync -av ${TMPDIR}/out/* ${pathToZeno}/CLARK_out/virulence/31kmer

echo -e "\nFinished running clark01.sh at `date`."
