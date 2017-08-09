#!/bin/bash
#$ -cwd
#$ -j n
#$ -l h_rt=96:00:00
#$ -m bea
#$ -M harryjreed1985@gmail.com
#$ -N repair_PE7
#$ -pe smp 4
#$ -l os=rhel6.3
#$ -l vf=2G
#$ -l h_vmem=4G
#$ -l zenodotus=true

# Path variables
export pathFromZeno=/zenodotus/masonlab/rotcollab_scratch/nbo5/trim_amb_samples/trim_amb_adapt
export pathToZeno=/zenodotus/masonlab/rotcollab_scratch/har2011/ambulance/FASTQ_hudson_alpha/repair_PEs
export pathToSoft=/home/har2011/bin

# Sync data to tmpdir
echo -e "\nrepair_PE7.sh started on ambulance samples at `date`.\nSyncing from: ${pathZeno} \nTo: ${TMPDIR}/working\n"
cd $TMPDIR
mkdir completed
mkdir working
mkdir soft

cd ${pathFromZeno}
rsync -av $( ls | egrep "AE|AS|AW" | sort | grep -o "A.*" | head -n 700 | tail -n 100 )  ${TMPDIR}/working
rsync -av ${pathToSoft}/parallel ${TMPDIR}/soft
rsync -av ${pathToSoft}/fastqCombinePairedEnd.py ${TMPDIR}/soft

# Unzipping files
echo -e "\nStarted to gunzip files at `date`.\n"
cd $TMPDIR/working
ls | time ${TMPDIR}/soft/parallel --gnu -j $NSLOTS gunzip

# renaming files
for reN in $( find . -name '*.clean' | sort | grep -o "A....."); do
        mv ${reN}_pair01.fastq.out.clean ${reN}_pair01.fastq
        mv ${reN}_pair02.fastq.out.clean ${reN}_pair02.fastq
done

# Create reference files for parallel
find . -name '*_pair01.fastq' | sort | grep -o 'A.*' > ref_pair01
find . -name '*_pair02.fastq' | sort | grep -o 'A.*' > ref_pair02
cat ref_pair01 ref_pair02

# run fastqCombinePairedEnd.py in parallel 
time ${TMPDIR}/soft/parallel --gnu --xapply -j $NSLOTS python ${TMPDIR}/soft/fastqCombinePairedEnd.py :::: ref_pair01 :::: ref_pair02

# Mv all completed files to completed
for mvRen in $( find . -name '*.fastq' | sort | grep -o "A....."); do
	mv ${mvRen}_pair01.fastq_pairs_R1.fastq ../completed/${mvRen}_pairs_R1.fastq
	mv ${mvRen}_pair02.fastq_pairs_R2.fastq ../completed/${mvRen}_pairs_R2.fastq
	mv ${mvRen}_pair01.fastq_singles.fastq ../completed/${mvRen}_singles.fastq
done

# Gnzip files and sync back to zeno
cd ${TMPDIR}/completed
find . -name '*.fastq' | sort | grep -o 'A.*' > gzRef
time ${TMPDIR}/soft/parallel --gnu -j $NSLOTS gzip :::: gzRef
rsync -av $( ls | egrep "*.gz" | sort | grep -o "A.*" )  ${pathToZeno}

echo -e "\nFinished repair_PE7.sh at `date`."
