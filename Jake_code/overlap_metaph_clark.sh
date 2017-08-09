#!/bin/bash
# overlap_metaph_clark.sh
# Looks for overlaping species between clark abundance output and metaphlan. 

# Declaring variables
in_files=''
out_file=''
out_d=''

help_readout="\nUsage:overlap_metaph_clark.sh -i <input_files> -d <out_dir> -f <out_file>
	-h	--help		Displays this message
	-i	--input		Input files, separated by a comma
	-d	--output_dir	Output directory
	-f	--output_file	Output filename
	*	Long form of options require \"=\"
		"

# Check for file
if [[ -z $@ ]]; then
        echo -e "$help_readout"
        exit 1
fi

# GetOpt
ARGS=$( getopt -o h::i:d:f: -l "help::,input:,output_dir:,output_file:" -n "overlap_metaph_clark.sh" -- "$@" );

eval set -- "$ARGS";

# extract options and their arguments into variables
while true; do
        case "$1" in
                -h|--help)
                shift;
                        echo -e "${help_readout}";
                exit 1;
                ;;
                -i|--input)
                shift;
                if [[ -n $1 ]]; then
			in_files=$1;
		shift;
		fi
		;;
		-d|--output_dir)
		shift;
		if [[ -n $1 ]]; then
			out_d=$1;
		shift;
		fi
		;;
		-f|--output_file)
		shift;
		if [[ -n $1 ]]; then
			out_file=$1;
		shift;
		fi
		;;
		--)
		shift;
		break;
		;;
	esac
done

# Parsing input files 
in_clark=$( echo ${in_files} | cut -d ',' -f1 );
in_meta=$( echo ${in_files} | cut -d ',' -f2 );

echo -e "\nVariable echo:
	Input files: ${in_files}
	Input clark file: ${in_clark}
	Input metaphlan file: ${in_meta}
	Ouput directory: ${out_d}
	Ouput file name: ${out_file}
	"
# Search for overlapping species
echo "Species,Clark_rel_abun,Meta_rel_abun" > ${out_d}/${out_file}_spp
echo "Genus" > ${out_d}/${out_file}_gen
echo "tmp" > ${out_d}/${out_file}_gen.tmp
echo -e "\nProcessing files ${in_clark} and ${in_meta}."
while read -r line || [[ -n "$line" ]]; do
	clark_spp_nm=$( echo "$line" | cut -d ',' -f1 );
	if [[ $clark_spp_nm != "Name" && $clark_spp_nm != "UNKNOWN" ]]; then
		clark_spp_nm_meta_fmt=$( echo "$clark_spp_nm" | sed -e "s: :_:g" );
		overlap_chk_spp=$( egrep "s__${clark_spp_nm_meta_fmt}" ${in_meta} | head -n 1 );
		clark_gen_nm_meta_fmt=$( echo "$clark_spp_nm" | cut -d ' ' -f1 );
		overlap_chk_gen=$( egrep "g__${clark_gen_nm_meta_fmt}" ${in_meta} );
		
		if [[ -n $overlap_chk_spp ]]; then
			echo -e "\n\t\tFound species ${clark_spp_nm} in both files"
			clark_rel_abun=$( echo "$line" | cut -d ',' -f6 );
			meta_rel_abun=$( echo "$overlap_chk_spp" | cut -f2 );
			echo "${clark_spp_nm},${clark_rel_abun},${meta_rel_abun}" >> ${out_d}/${out_file}_spp
		fi
		if [[ -n $overlap_chk_gen ]]; then
			echo -e "\n\t\tFound genus ${clark_gen_nm_meta_fmt} in both files"
			echo "${clark_gen_nm_meta_fmt}" >> ${out_d}/${out_file}_gen.tmp
		fi
	fi		
done < ${in_clark}

tail -n +2 ${out_d}/${out_file}_gen.tmp | sort | uniq >> ${out_d}/${out_file}_gen
rm ${out_d}/${out_file}_gen.tmp

# End script
