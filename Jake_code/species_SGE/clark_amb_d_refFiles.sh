#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_rt=96:00:00
#$ -m bea
#$ -M harryjreed1985@gmail.com
#$ -N clark_dRef
#$ -pe smp 8
#$ -l os=rhel6.3
#$ -l vf=25G
#$ -l h_vmem=35G
#$ -l zenodotus=true

# Path variables
pathFromZeno=/zenodotus/masonlab/rotcollab_scratch/har2011/ambulance/FASTQ_hudson_alpha/repair_PEs/pairs
pathToZeno=/zenodotus/masonlab/rotcollab_scratch/har2011/ambulance
pathToSoft=/home/har2011/bin
pathToDB=/home/har2011/bin/clark_db

# Sync data to tmpdir
cd $TMPDIR

echo -e "\nclark_ambulance_debug_ref.sh started on ambulance samples at `date`.\nSyncing from: ${pathZeno} \nTo: ${TMPDIR}/clark01_work\n"

mkdir out
mkdir work
mkdir abundance_out
mkdir soft
mkdir soft/clark_db

# Testing on only largest files
rsync -av $( find ${pathFromZeno}/*.gz -size +2300M )  ${TMPDIR}/work
rsync -a ${pathToSoft}/CLARKSCV1.2.2-b ${TMPDIR}/soft
rsync -av ${pathToSoft}/parallel ${TMPDIR}/soft

# Sync only necessary DB files
rsync -a ${pathToDB} ${TMPDIR}/soft


# Unzipping files
echo -e "\nStarted to gunzip files at `date`.\n"
cd $TMPDIR/work
ls | ${TMPDIR}/soft/parallel --gnu -j $NSLOTS gunzip
echo -e "Finished gunzip files at `date`.\n"

# Setting targets from database
cd $TMPDIR/soft/CLARKSCV1.2.2-b 
./set_targets.sh ${TMPDIR}/soft/clark_db bacteria --species

# Create ref files for clark
cd $TMPDIR/work
find ./*pairs_R1.fastq | sort > refPair01.txt
find ./*pairs_R2.fastq | sort > refPair02.txt

find ./*pairs_R1.fastq | sort > refResults.txt
sed -i 's:./:${TMPDIR}/out/:g' refResults.txt
sed -i 's:pairs_R1.fastq.clean:_clarkOut:g' refResults.txt
echo -e "Following is cat of refResults.txt:\n"
cat refResults.txt

# Classifying metagenome PE reads
echo -e "\nMetagenome classification began at `date`.\n"
cd $TMPDIR/soft/CLARKSCV1.2.2-b 	
./classify_metagenome.sh -P ${TMPDIR}/work/refPair01.txt ${TMPDIR}/work/refPair02.txt -R ${TMPDIR}/work/refResults.txt -m 0 -n $NSLOTS --ldm -k 25
echo -e "\nMetagenome classification completed at `date`.\n"

# Abundance estimation of output files with high-confidence
echo -e "Abundance estimation began at `date`.\n"
cd $TMPDIR/soft
for clarkFile in $( find ${TMPDIR}/out/*.csv )
do
	./estimate_abundance.sh -F ${clarkFile} -D ${pathToDB} --highconfidence -a 0.001 > ${clarkFile}_abundance
done

# Mv files and rename to clark_01_abundance_out
mv ${TMPDIR}/out/*_abundance ${TMPDIR}/abundance_out/
cd ${TMPDIR}/abundance_out/
for rename in $( ls | egrep "AE|AS|AW" | grep -o 'A.....' ) 
do 
	mv ${rename}_clarkOut.csv_abundance ./${rename}_abundance.csv 
	echo -e "\nRenaming file ${rename}."
done
echo -e "\nAbundance estimation and renaming finished at `date`.\n"

# Zipping clark01_out files
echo -e "Zipping files started at `date`.\n"
cd $TMPDIR/out
find ./*.csv | sort > gzRef 
${TMPDIR}/soft/parallel --gnu -j $NSLOTS gzip :::: gzRef
rm gzRef
echo -e "\nFinished zipping files at `date`.\n"

# Sync data back to $myamdata
echo -e "Syncing files back to zeno.\n"
rsync -av ${TMPDIR}/out/* ${pathToZeno}/CLARK_out
rsync -av ${TMPDIR}/abundance_out/* ${pathToZeno}/CLARK_abundance_out

echo -e "\nFinished running clark_ambulance_debug_ref.sh at `date`."
