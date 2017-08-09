#!/bin/bash
# JR_estimate_abundance.sh
# Estimate abundance from custom database 

# Starting variables
in_sample=''
out_f=''
targets_f=''
out_d=''
filter_score=''

readout_inf="
Usage:JR_estimate_abundance.sh -i <input_sample> -t <targets_file> -o <output_file> -d <out_dir> -s <filter_score>

	-h|--help		Displays this message
	-i|--input		Input clark sample output
	-t|--targets		Targets file used for clark classification
	-o|--outFile		Output file
	-d|--outDir		Output directory
	-s|--filterScore	Filter Score, should be between 1 and 100, def is 30
	
	"

# Test for arguments
if [[ -z $1 ]]; then
	echo -e "${readout_inf}"
	exit 1
fi

#Set options
ARGS=$( getopt -o h::i:t:o:d:s: -l "help::,input:,targets:,outFile:,outDir:,filterScore" -n "JR_estimate_abundance.sh" -- "$@" );

eval set -- "$ARGS";
while true; do
        case "$1" in
                -h|--help)
                        shift;
                        echo -e "${readout_inf}";
                        exit 1;
                ;;
                -i|--input)
                        shift;
                        if [[ -n $1 ]]; then
                                in_sample=$1;
                                echo ${in_sample}
                        shift;
                        fi
                ;;
 		-t|--targets)
                        shift;
                        if [[ -n $1 ]]; then
                                targets_f=$1;
                                echo ${targets_f}
                        shift;
                        fi
                ;;
		-o|--outFile)
                        shift;
                        if [[ -n $1 ]]; then
                                out_f=$1;
                                echo ${out_f}
                        shift;
                        fi
                ;;
		-d|--outDir)
			shift;
			if [[ -n $1 ]]; then
				out_d=$1;
				echo ${out_d}
			shift;
			fi
		;;
		-s|--filterScore)
			shift; 
			if [[ -n $1 ]]; then
				filter_score=$1;
				echo ${filter_score}
			shift;
			fi
		;;
         	--)
                        shift;
                        break;
                ;;
        esac
done

in_sample_nm=$( echo -e "${in_sample}" | sed -e "s:.*/::;s/$//"  );
out_f_nm=$( echo -e "${out_f}" | sed -e "s:.*/::;s/$//"  );

echo -e "
Input clark out file: ${in_sample}
  name: ${in_sample_nm}
Output file and path: ${out_f}
  name: ${out_f_nm}
Targets file: ${targets_f}
Output directory: ${out_d}
Filter score: ${filter_score}
	"

# Checks
if [[ -z $in_sample || -z $out_f || -z $targets_f || -z $out_d ]]; then
	echo -e "
Required arguments were not met; -i,-f,-t,-d must be supplied, exiting
		"
	exit 1;
fi

if [[ -f $in_sample && -d $out_d && -f $targets_f ]]; then
	echo -e "
Required inputs meet file/directory standards, proceeding.
		"
else 
	echo -e "
Required inputs do not meet file/directory standards, exiting.
		"
	exit 1;
fi

# Begin parsing
tot_reads=$( tail -n +2 ${in_sample} | wc -l );
tot_antibiotic_markers=$( cat ${targets_f} | wc -l );

# Sort and only keep files without NA for 1st assignment score > 30, fairly stringent
echo -e "
Starting to sort input clark classification output based on input parameters at `date`
	"

if [[ -z $filter_score ]]; then
	cat ${in_sample} | awk '$4 != "NA"' FS=',' | awk '$5 > 30' FS=',' > ${out_d}/${in_sample_nm}.tmpfile_score
else
	cat ${in_sample} | awk '$4 != "NA"' FS=',' | awk -v var="${filter_score}" '$5 > var' FS=',' > ${out_d}/${in_sample_nm}.tmpfile_score
fi

# for loop for processing all virulence markers, relative abundance output in percentage as in 95.0%
echo -e "
Starting to process each assignment based on targets and sample file based on input parameters at `date`
	"
echo -e "Antibiotic_resistance_marker,Number_hits,Relative_abundance" > ${out_d}/${out_f_nm}.tmp

for vir in $( seq 1 ${tot_antibiotic_markers} ); do
	num_hits=$( cat ${out_d}/${in_sample_nm}.tmpfile_score | awk -v var=${vir} '$4 == var' FS=',' | wc -l );
	rel_abun=$( echo "scale=10; (${num_hits}/${tot_reads})*100" | bc );
	vir_marker=$( awk -v var=${vir} '$2 == var {print $1}' FS='\t' ${targets_f} ); 
	vir_marker_nm=$( echo "${vir_marker}" | sed -e "s:.*/::;s/$//" | sed -e "s:.fa::g" );
	echo -e "${vir_marker_nm},${num_hits},${rel_abun}" >> ${out_d}/${out_f_nm}.tmp
done

# Creating final file
head -n 1 ${out_d}/${out_f_nm}.tmp > ${out_d}/${out_f_nm}
tail -n +2 ${out_d}/${out_f_nm}.tmp | awk '$2 != 0' FS=',' | sort -nrk2 -t ','  >> ${out_d}/${out_f_nm}

rm ${out_d}/${out_f_nm}.tmp
rm ${out_d}/${in_sample_nm}.tmpfile_score

echo -e "
Finished processing abundance at `date`
	" 
