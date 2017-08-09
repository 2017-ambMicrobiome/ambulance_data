#!/bin/bash
# recat 1st 100 files for ambulance project

# navigate to dir
cd /zenodotus/masonlab/pathomap_scratch/ebrahim/Ambulance_Study/hudson_alpha/FASTQ

# create cross reference file 
dir | egrep "AE|AS|AW" | head -n 100 | egrep "_01" | grep -o "A....." > /home/har2011/SGE_scripts/catPairedEnds/recat100/recat100CrossRef

# navigate to scrip dir
cd /home/har2011/SGE_scripts/catPairedEnds

# create ref file from catPair01 and catPair02 
touch ./recat100/RefFile100
while read line; do
cat catPair01_ambulance.txt | egrep "${line}" >> ./recat100/RefFile100
cat catPair02_ambulance.txt | egrep "${line}" >> ./recat100/RefFile100
done < ./recat100/recat100CrossRef

# End recat100.sh
