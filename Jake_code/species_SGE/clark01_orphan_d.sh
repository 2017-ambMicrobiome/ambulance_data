#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_rt=96:00:00
#$ -m bea
#$ -M harryjreed1985@gmail.com
#$ -N clark01_orphan_d
#$ -pe smp 8
#$ -l os=rhel6.3
#$ -l vf=20G
#$ -l h_vmem=300G
#$ -l zenodotus=true

# Path variables
pathFromZeno=/zenodotus/masonlab/rotcollab_scratch/har2011/ambulance/FASTQ_hudson_alpha/repair_PEs/orphans
pathToZeno=/zenodotus/masonlab/rotcollab_scratch/har2011/ambulance
pathToSoft=/home/har2011/bin
pathToDB=/home/har2011/bin/clark_db

# Sync data to tmpdir
cd $TMPDIR

echo -e "\nclark01_orphan.sh started on ambulance samples at `date`.\nSyncing from: ${pathZeno} \nTo: ${TMPDIR}/clark01_work\n"

mkdir out
mkdir work
mkdir abundance_out
mkdir soft
mkdir soft/clark_db

# 0-200
rsync -av $( find ${pathFromZeno}/*.gz | sort | head -n 2 )  ${TMPDIR}/work
rsync -a ${pathToSoft}/CLARKSCV1.2.2-b ${TMPDIR}/soft
rsync -av ${pathToSoft}/parallel ${TMPDIR}/soft

# Sync only necessary DB files
rsync -a ${pathToDB} ${TMPDIR}/soft

# Unzipping files
echo -e "\nStarted to gunzip files at `date`.\n"
cd $TMPDIR/work
ls | ${TMPDIR}/soft/parallel --gnu -j $NSLOTS gunzip
echo -e "Finished gunzip files at `date`.\n"

# Making Ref files 
cd $TMPDIR/work
find `pwd` -name '*.fastq' | sort > singlesRef_in.txt
cp singlesRef_in.txt singlesRef_out.txt
sed -i "s:/work/:/out/:g" singlesRef_out.txt
sed -i "s:.fastq:_clarkOut:g" singlesRef_out.txt
cat singlesRef_in.txt
cat singlesRef_out.txt
ls

# Setting targets from database
cd $TMPDIR/soft/CLARKSCV1.2.2-b 
./set_targets.sh ${TMPDIR}/soft/clark_db bacteria --phylum

# Classifying metagenome Singles
cd $TMPDIR/soft/CLARKSCV1.2.2-b
./classify_metagenome.sh -O ${TMPDIR}/work/singlesRef_in.txt -R ${TMPDIR}/work/singlesRef_out.txt -m 0 -n $NSLOTS --ldm -k 25
rm $TMPDIR/work/singlesRef_in.txt
rm $TMPDIR/work/singlesRef_out.txt

# Determine format of clark output file
echo -e 'Creating formatClarkOut to look at output file types'
cd $TMPDIR/out
find . > formatClarkOut

# Abundance estimation of output files with high-confidence
echo -e "Abundance estimation began at `date`.\n"
cd $TMPDIR/soft/CLARKSCV1.2.2-b
for clarkFile in $( find ${TMPDIR}/out/*.csv )
do
	./estimate_abundance.sh -F ${clarkFile} -D ${pathToDB} -c 0.95 -g 0.06 -a 0.001 > ${clarkFile}_abundance
done

# Mv files and rename to clark_01_abundance_out
mv ${TMPDIR}/out/*_abundance ${TMPDIR}/abundance_out/
cd ${TMPDIR}/abundance_out/
for rename in $( ls | egrep "AE|AS|AW" | grep -o 'A.....' ) 
do 
	mv ${rename}_singles_clarkOut.csv_abundance ./${rename}_singles_abundance.csv 
	echo -e "\nRenaming file ${rename}."
done
echo -e "\nAbundance estimation and renaming finished at `date`.\n"

# Zipping clark01_orphan_out files
echo -e "Zipping files started at `date`.\n"
cd $TMPDIR/out
find ./*.csv | sort > gzRef 
${TMPDIR}/soft/parallel --gnu -j $NSLOTS gzip :::: gzRef
rm gzRef
echo -e "\nFinished zipping files at `date`.\n"

# Sync data back to $myamdata
echo -e "Syncing files back to zeno.\n"
rsync -av ${TMPDIR}/out/* ${pathToZeno}/CLARK_out_orphans
rsync -av ${TMPDIR}/abundance_out/* ${pathToZeno}/CLARK_abundance_out/estimate_param_fromEval

echo -e "\nFinished running clark01_orphan.sh at `date`."
