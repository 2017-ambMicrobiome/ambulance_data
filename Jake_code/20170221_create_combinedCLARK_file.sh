#!/bin/bash
# create_combined file from clark data

# Dir variables
data_d=/zenodotus/masonlab/rotcollab_scratch/har2011/metagenomics/CLARK_abundance_out/param_fromEval_species
tmp_d=/zenodotus/masonlab/rotcollab_scratch/har2011/tmp

# Sync only pair data
rsync -av $(find ${data_d} | egrep "*pairs*") ${tmp_d}

# Loop over files and create tsv documents for processing with merge_metaphlan
cd ${tmp_d}
for csv in ./*.csv; do
	echo -e "\nProcessing file ${csv} for use in merge_metaphlan at `date`"
	tail -n +2 ${csv} | grep -v "UNKNOWN" | cut -d',' -f1,5 | sed 's/,/\t/g' > ${csv}.tsv
	egrep "UNKNOWN" ${csv} | cut -d',' -f1,4 | sed 's/,/\t/g' >> ${csv}.tsv
done

# Use metaphlan merge to create output file
python /home/har2011/software/metaphlan/utils/merge_metaphlan_tables.py $(find . -name "*.tsv") > 20170221_AllClark_combined.tsv

# End Script
