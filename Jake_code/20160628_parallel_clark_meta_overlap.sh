#!/bin/bash
# 20160628_parallel_clark_meta_overlap.sh

export tmp=$TMP
cd $tmp
mkdir raw_f
mkdir ref
mkdir out
export out=${tmp}/out
export raw_f=${tmp}/raw_f
export ref=${tmp}/ref

parallel -j 24 --gnu rsync -av {1} ${raw_f} ::: $( find /zenodotus/masonlab/pathomap_scratch/ebrahim/Ambulance_Study/clean_up/output/*.out )
parallel -j 24 --gnu rsync -av {1} ${raw_f} ::: $( find /zenodotus/masonlab/rotcollab_scratch/har2011/ambulance/CLARK_abundance_out/estimate_param_fromEval/*pairs* );
sleep 2s

touch ${ref}/allClark.fofn
touch ${ref}/allMeta.fofn
find ${raw_f}/*.out >> ${ref}/allMeta.fofn
find ${raw_f}/*_pairs_abundance.csv >> ${ref}/allClark.fofn

end_seq=$( cat ${ref}/allClark.fofn | wc -l );
sleep 2s

touch ${ref}/input.fofn
for file in $( seq 1 $end_seq ); do
	clark_f=$( sed -n "${file}p" ${ref}/allClark.fofn );
	meta_f=$( sed -n "${file}p" ${ref}/allMeta.fofn );
	echo "${clark_f},${meta_f}" >> ${ref}/input.fofn
done

parallel -k -j 24 --gnu overlap_metaph_clark.sh -i {} -d ${out} -f {/}_overlap.csv :::: ${ref}/input.fofn

rsync -av ${out}/* /zenodotus/masonlab/rotcollab_scratch/har2011/ambulance/clark_meta_overlap	

cd $tmp
rm -r *
