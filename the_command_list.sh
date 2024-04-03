#
#mkdir basic_info/
vcftools --gzvcf $2/$1.vcf.gz --remove-filtered-all --hwe 1e-6 --max-missing 1 --counts2 --derived --max-alleles 2 --min-alleles 2--remove-indels --out basic_info/$1 &&  ## --derived must work 
vcftools --gzvcf $2/$1.vcf.gz --remove-filtered-all --hwe 1e-6 --max-missing 1 --get-INFO AA --max-alleles 2 --min-alleles 2 --remove-indels --out basic_info/$1 &&
python /media/sheldon/huchy/mosuo/mosuo_latest_2020_10_22/scripts/filteraa_calls.py basic_info/ $1 && # get $1.frq.count.uppercaseaa && $1.INFO.uppercaseaa
less basic_info/$1.INFO.uppercaseaa | cut -f1,2 | sed 1d > basic_info/$1.pos.uppercaseaa &&
vcftools --gzvcf $2/$1.vcf.gz --positions basic_info/$1.pos.uppercaseaa --extract-FORMAT-info GT --remove-indels --out basic_info/$1.uppercaseaa &&
python /media/sheldon/huchy/mosuo/mosuo_latest_2020_10_22/scripts/ac_find_indivs_data.py basic_info $1 && # get $1.alt.count.indiv
python /media/sheldon/huchy/mosuo/mosuo_latest_2020_10_22/scripts/write_dots.py basic_info/$1.INFO.uppercaseaa &&
mkdir $1 &&
python /media/sheldon/huchy/mosuo/mosuo_latest_2020_10_22/scripts/mac_find_indivs_data.py basic_info/ $1 uppercaseaa &&
singularity exec -B /media/sheldon/huchy/program/.vep:/home/.vep,/media/sheldon/huchy/ani_from_zhou/20240227/6mammals/basic_info:/home/basic_info vep.sif vep --cache --offline --species $3 --dir_cache /home/.vep --format vcf -i /home/basic_info/$1.INFO.uppercaseaa.vepinput -o /home/basic_info/$1.vep_output.txt --everything --pick --force_overwrite --fork 60 --cache_version ${5:-101} &&
## --cache_version 101 by default
## --cache_version 97 for Mmul_8.0.1
## --cache_version 90 for gorGor3.1
## --cache_version 101 for gorGor4
###sed 's/\t/:/g' basic_info/$1.pos.uppercaseaa >basic_info/$1.pos.uppercaseaa_1 &&
###awk 'NR==FNR{a[$1];next}$2 in a' basic_info/$1.pos.uppercaseaa_1 vcf_files/$3.vep_output.txt >$1/$1.vep_output.txt &&
grep -v "#" basic_info/$1.vep_output.txt | cut -f 2,4,7 > $1/$1.functype.vep.txt &&
cp basic_info/$1.frq.count.indiv.uppercaseaa $1/ &&
python /media/sheldon/huchy/mosuo/mosuo_latest_2020_10_22/scripts/merge_variants_types_given_vep.py $1/ $1 vep &&
gzip basic_info/$1.uppercaseaa.GT.FORMAT &&
cd $1/ &&
{ 
cat $1.frq.count.indiv.functype.vep.uppercaseaa | grep 'missense_variant' > $1.frq.count.indiv.missense.vep.uppercaseaa
cat $1.frq.count.indiv.functype.vep.uppercaseaa | grep 'stop_gained' > $1.frq.count.indiv.stopgain.vep.uppercaseaa
cat $1.frq.count.indiv.functype.vep.uppercaseaa | grep 'stop_lost' > $1.frq.count.indiv.stoplost.vep.uppercaseaa
cat $1.frq.count.indiv.stoplost.vep.uppercaseaa $1.frq.count.indiv.stopgain.vep.uppercaseaa > $1.frq.count.indiv.nonsense.vep.uppercaseaa
cat $1.frq.count.indiv.functype.vep.uppercaseaa | grep 'synonymous_variant' > $1.frq.count.indiv.codingsynon.vep.uppercaseaa 
cat $1.frq.count.indiv.functype.vep.uppercaseaa | grep "synonymous_variant" | grep -v "synonymous_variant,"| grep -v ",synonymous_variant" > $1.frq.count.indiv.codingsynononly.vep.uppercaseaa
cat $1.frq.count.indiv.functype.vep.uppercaseaa | grep -w "splice_donor_variant" | grep -v "splice_donor_variant," > $1.frq.count.indiv.splicedonor.vep.uppercaseaa
cat $1.frq.count.indiv.functype.vep.uppercaseaa | grep -w "splice_acceptor_variant" | grep -v "splice_acceptor_variant," > $1.frq.count.indiv.spliceacceptor.vep.uppercaseaa
cat $1.frq.count.indiv.splicedonor.vep.uppercaseaa $1.frq.count.indiv.spliceacceptor.vep.uppercaseaa > $1.frq.count.indiv.splice.vep.uppercaseaa
cat $1.frq.count.indiv.splice.vep.uppercaseaa $1.frq.count.indiv.nonsense.vep.uppercaseaa > $1.frq.count.indiv.splicenonsense.vep.uppercaseaa
} &&
cd ../ &&
for ii in {"splicenonsense","missense","codingsynononly"}
do
	#python /media/sheldon/huchy/ani_from_zhou/slim/2.mac_count_all_data_byloci_chr.py $4 $1/ $1 vep $ii $1/
	python /media/sheldon/huchy/ani_from_zhou/slim/2.mac_count_all_data_byloci_chr.py samples/$4.sample.txt $1/ $1 vep $ii $1/ 
	#gzip $1/$1.frq.count.indiv.$ii.vep.uppercaseaa 
done
echo done


