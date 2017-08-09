#!/bin/bash
# 20160822_convert_meta_clark_overlap_spp_biom.sh
# Used for standardizing meta by total reads and converting to biom format for vegan package processing in R

mkdir ${mydata_main}/tmp/work

# Dir variables
export conv=/home/har2011/software/metaphlan/conversion_scripts
export utils=/home/har2011/software/metaphlan/utils
export work=${mydata_main}/tmp/work
export bk=${amb}/clark_meta_overlap/processed_files/spp
export ref=${amb}/referenceFiles
export tmp=${mydata_main}/tmp
export out=${amb}/diversity/clark_meta_overlap_spp

# Copying bk files
rsync -av $( find $bk/*_spp ) $work

# Replacing all headers with file name
cd $work
for file in $( find ./*_spp ); do
	
	chk=$( cat ${file} )
	if [[ -z ${chk} ]]
	then
		echo -e "\nMoving ${file} up a dir because it was empty."
		mv ${file} ../
	fi
	
	echo -e "Renaming header for file ${file}\n"
	fil_name=$( echo ${file} | grep -o 'A.....' )
	sed -i "s:#SampleID:ID:g" $file
	sed -i "s:Metaphlan2_Analysis:${fil_name}:g" $file
	tot_reads=$( egrep "${fil_name}" ${ref}/ambulanceMetaComplete.csv | cut -d',' -f10 )
	
	if [[ $tot_reads != "NA" ]]
	then
		echo -e "The total reads for ${fil_name} are ${tot_reads}\n" 
		python ${conv}/metaphlan2biom.py --nreads ${tot_reads} ${file} > ${fil_name}.biom

		echo -e "\nConverting ${fil_name}.biom to *.tsv"
	
		biom convert -i ${fil_name}.biom -o ${fil_name}.tsv --to-tsv --header-key taxonomy
		awk -v OFS='\t' -F'\t' '{print $3, $2}' ${fil_name}.tsv > ${fil_name}.meta.out
		sed -i "s:; :|:g" ${fil_name}.meta.out
		sed -i '1,2d' ${fil_name}.meta.out
		mv ${fil_name}.meta.out ${fil_name}
	else
		echo -e "\n moving ${file} up a dir due to total reads equaling 'NA'"
		mv ${file} ../
	fi
done

# Cleaning up
rm *.out *.biom *.tsv

# Combining and starting final conversions
echo -e "Merge metaphlan tables started at `date`\n"
python ${utils}/merge_metaphlan_tables.py * > amb_clark_meta_ovr_spp.out

# Merge metaphlan tables converting to biom format
echo -e "Merged metaphlan tables converting to biom started at `date`\n"
python ${conv}/metaphlan2biom.py --nreads 1 amb_clark_meta_ovr_spp.out > amb_clark_meta_ovr_spp.biom

# At this point all counts are num_read hits/100
echo -e "Final conversion to tsv document started at `date`\n"
biom convert -i amb_clark_meta_ovr_spp.biom -o amb_clark_meta_ovr_spp.biom.tsv --to-tsv --header-key taxonomy
mv amb_clark* ../

# clean up 
rm *

# Sed to remove headers from tsv to use for vegan diversity in R 
cd ${tmp}
sed -i '1d' amb_meta.biom.tsv
sed -i 's:#::g' amb_meta.biom.tsv
rev amb_clark_meta_ovr_spp.biom.tsv | cut -f 2- | rev > amb_clark_meta_ovr_spp.biom.tsv_2
rm amb_clark_meta_ovr_spp.biom.tsv; mv amb_clark_meta_ovr_spp.biom.tsv_2 amb_clark_meta_ovr_spp.biom.tsv
mv amb* ${out}
