#!/bin/bash
# 20160805_parallel_custom_abundance_ebr.sh

# Create directories

export tmp=/scratchLocal/cmlab/har2011/tmp_2
cd $tmp
mkdir raw_f
mkdir ref_f
mkdir abun
export ref_f=${tmp}/ref_f
export raw_f=${tmp}/raw_f
export abun_o=${tmp}/abun
export script_d=/home/har2011/software/CLARKSCV1.2.2-b

# Sync files generated from clark using custom CARD database
parallel -j 24 --gnu rsync -av {1} ${raw_f} ::: $( find /zenodotus/masonlab/rotcollab_scratch/har2011/metagenomics/CLARK_out/clark_card_ambulance/31kmer/*.gz )
parallel -j 24 --gnu gunzip {1} ::: $( find ${raw_f}/*.gz )

# Create file reference
touch ${ref_f}/clarkOut.fofn
find ${raw_f}/*.csv >> ${ref_f}/clarkOut.fofn

# Process files
parallel -k -j 24 --gnu ${script_d}/JR_estimate_abundance_all.sh -i {}  -d ${abun_o} -o {/}_abundance.csv -s 50 -t /zenodotus/masonlab/rotcollab_scratch/har2011/DB_files/clark_virulence_db/targets_finished.txt :::: ${ref_f}/clarkOut.fofn

# Rsync files back to abundance directory
parallel -j 24 --gnu rsync -av {1} /zenodotus/masonlab/rotcollab_scratch/har2011/metagenomics/CLARK_abundance_out/clark_card_ambulance/31kmer/score_50_all ::: $( find ${abun_o}/*.csv )
rm -r ${tmp}/*
