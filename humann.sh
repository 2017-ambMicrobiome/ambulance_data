#$ -e /zenodotus/masonlab/pathomap_scratch/ebrahim/Ambulance_Study/humann_output/AW0508_clean.8.2.humann2.err
#$ -m bea
#$ -M eba2001@med.cornell.edu
#$ -l zenodotus=true
#$ -l h_vmem=2G
#$ -pe smp 8
#$ -N AW0508_clean.8.2.humann2

name=AW0508_clean.8.2
out_name=AW0508_clean.humann2
input_file=/zenodotus/masonlab/pathomap_scratch/ebrahim/Ambulance_Study/clean_up/AW0508_clean.fastq.gz
cpu=8
out_dir=/zenodotus/masonlab/pathomap_scratch/ebrahim/Ambulance_Study/humann_output

rsync -av $input_file $TMPDIR/
cd $TMPDIR

export PATH=/home/emh2013/programs/jdk1.8.0_66/bin:/home/emh2013/programs/jdk1.8.0_66:/home/emh2013/anaconda/bin:/home/emh2013/programs/samtools-1.3/bin:/usr/lib64/qt-3.3/bin:/usr/kerberos/sbin:/usr/kerberos/bin:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/opt/bin:/opt/dell/srvadmin/bin:/home/ebrahim/bin:/home/yos2006/tools/BEDTools-Version-2.14.3/bin:/home/emh2013/scripts/bioinfo_scripts/python:/home/emh2013/scripts/bioinfo_scripts/perl:/home/emh2013/scripts/bioinfo_scripts/awk:/home/emh2013/scripts/bioinfo_scripts/bash:/home/emh2013/scripts/bioinfo_scripts/scripts_CRG:/home/emh2013/dev/jitterbug-code:/home/emh2013/dev/jitterbug-code/jip_scripts/:/home/emh2013/dev/jitterbug-code/scripts/:/home/emh2013/programs/htop-1.0/bin/bin:/home/emh2013/programs/FastQC:/home/emh2013/programs/lumpy-sv/scripts:/home/emh2013/programs/lumpy-sv/bin:/home/emh2013/programs/IGVTools:/home/emh2013/programs/breakdancer-1.1_2011_02_21/cpp:/home/emh2013/programs/scalpel-0.3.1:/home/emh2013/programs/ncbi-blast-2.2.30+/bin:/home/emh2013/programs/last-572/src:/home/emh2013/programs/last-572/scripts:/home/emh2013/programs/nanopore-scripts:/home/emh2013/programs/cdhit:/home/emh2013/programs/biobakery-shortbred-ddce9103c5ee:/home/emh2013/programs:/home/emh2013/programs/graphlan:/home/emh2013/programs/metaphlan2/utils:/home/emh2013/programs/metaphlan2/utils/export2graphlan:/home/emh2013/programs/metaphlan2:/home/emh2013/programs/lefse:/home/emh2013/programs/bowtie2-2.2.5:/home/emh2013/programs/humann2_v0.2.0/humann2/tools:/home/emh2013/programs/htslib-1.3/bin

echo "name: $name"
date=$(date)
echo " starting humann2 $name ... $date" >> /home/ebrahim/JOB_LOG.txt
humann2 --input $input_file  --threads $cpu  --remove-temp-output  --output-basename $out_name --output . 
date=$(date)
echo "DONE humann2 $name ... $date" >> /home/ebrahim/JOB_LOG.txt

echo "files generated:"
ls -l


rsync -av $out_name*   $out_dir/