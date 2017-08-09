#!/bin/bash
# mv files from /home/har2011 to /zenodotus/masonlab/rotcollab_scratch/har2011/FASTQ_hudson_alpha

mv -v /home/har2011/*.gz /zenodotus/masonlab/rotcollab_scratch/har2011/ambulance/FASTQ_hudson_alpha/

# Navigating to dir
cd /zenodotus/masonlab/rotcollab_scratch/har2011/ambulance/FASTQ_hudson_alpha/

# Create RefFile
ls | egrep "*.gz" > ./refMvFile.txt

# Create directories and mv files to new dir
for f in {1..8}
do 
	mkdir ./files${f}00; mv `cat ./refMvFile.txt | head -n ${f}00 | tail -n 100` ./files${f}00
done

# done mvAmbulanceFiles.sh
