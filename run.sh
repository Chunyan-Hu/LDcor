
## crossmap, bedtools
##*.NF.pVCF.genome.2nd_degree_unrelated : File contains the maximally sized set of unrelated individuals as determined from the relatedness analysis, - excludes all first and second degree relationships

## 工作路径及软链接
cd /path/to/LD_cor/
ln -s /path/to/EINSTEIN/ ./
ln -s /path/to/FL3_CL/ ./
ln -s /path/to/.vep/ ./

## 无亲源个体
less FL3_CL/FL3-CL-Dem.txt |grep "LGP Proband" |awk -F, '{print $2}' > LGP_Proband.list
less FL3_CL/FL3-CL-Dem.txt |grep "LGP Control Only" |awk -F, '{print $2}' > LGP_Control_Only.list
less EINSTEIN/relationships/EINSTEIN_Freeze_Three.NF.pVCF.genome.FILTERED_2nd-degree_relationships.genome |sed 1d |awk '{print$2}'|sort |uniq >NF.pVCF.genome.2nd_degree_unrelated.list_1
less EINSTEIN/relationships/EINSTEIN_Freeze_Three.NF.pVCF.genome.FILTERED_2nd-degree_relationships.genome |sed 1d |awk '{print$4}'|sort |uniq >NF.pVCF.genome.2nd_degree_unrelated.list_2
cat NF.pVCF.genome.2nd_degree_unrelated.list_1 NF.pVCF.genome.2nd_degree_unrelated.list_2|sort|uniq >NF.pVCF.genome.2nd_degree_unrelated.list
rm NF.pVCF.genome.2nd_degree_unrelated.list_1 NF.pVCF.genome.2nd_degree_unrelated.list_2
sort NF.pVCF.genome.2nd_degree_unrelated.list LGP_Proband.list |uniq -d > LGP_Proband.2nd_unrelated.list
sort NF.pVCF.genome.2nd_degree_unrelated.list LGP_Control_Only.list |uniq -d > LGP_Control_Only.2nd_unrelated.list

## 祖先基因型文件
#wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/ancestral_alignments/human_ancestor_GRCh37_e59.tar.bz2 ## AA.fasta
tar -jxvf human_ancestor_GRCh37_e59.tar.bz2
for ii in {1..22};do sed -i -e '1s/>ANCESTOR_for_chromosome:GRCh37:/>/g' -e '1s/:.*//g' human_ancestor_GRCh37_e59/human_ancestor_${ii}.fa ;done

## 坐标转换文件
#wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz ## chain file
less hg38ToHg19.over.chain.gz|sed 's/chr//g' |gzip >hg38ToHg19.over.chain.1.gz ## chr1 -> 1
echo '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele. Format: AA|REF|ALT|IndelType. AA: Ancestral allele, REF:Reference Allele, ALT:Alternate Allele, IndelType:Type of Indel (REF, ALT and IndelType are only defined for indels)">' > header.AA.txt ## header file

## 获取无亲源且插入祖先基因型的vcf, by chr
mkdir EINSTEIN_filter/ EINSTEIN_hg19/ EINSTEIN_AA/ EINSTEIN_withAA/ EINSTEIN_withAA/Proband/ EINSTEIN_withAA/Control/ output_Proband/ output_Control/
for ii in {1..22}
do
	{
	## make sure that the vcf file version is ok
	#bcftools view -v snps EINSTEIN/GL_by_chrom/EINSTEIN_Freeze_Three.${ii}.GL.vcf.gz -m2 -M2 |sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/g' |gzip -c > EINSTEIN_filter/chr${ii}.GL.hg38.snp.vcf.gz &&
	vcftools --gzvcf EINSTEIN/GL_by_chrom/EINSTEIN_Freeze_Three.${ii}.GL.vcf.gz --min-alleles 2 --max-alleles 2 --hwe 1e-15 --remove-indels --remove-filtered-all --recode --recode-INFO-all --stdout |gzip -c > EINSTEIN_filter/chr${ii}.GL.hg38.snp.vcf.gz && ## filter
	less -S EINSTEIN_filter/chr${ii}.GL.hg38.snp.vcf.gz |grep -v '^#'|awk -v OFS='\t' '{print $1,$2-1,$2}' > EINSTEIN_filter/chr${ii}.GL.hg38.snp.pos.bed && ## snps, (bed file is 0-based, [,))(5,9 in bed = [5,9) in bed = [6,10) in fa = (5,9] in fa)
	CrossMap.py bed hg38ToHg19.over.chain.1.gz EINSTEIN_filter/chr${ii}.GL.hg38.snp.pos.bed > EINSTEIN_hg19/chr${ii}.GL.hg38ToHg19.snp.pos && ## liftover 
	less -S EINSTEIN_hg19/chr${ii}.GL.hg38ToHg19.snp.pos |grep -v 'Unmap' |awk '{if($1==$5) print $5,$6,$7}' >EINSTEIN_hg19/chr${ii}.GL.hg19.snp.pos.bed && ## snp pos in hg19
	less -S EINSTEIN_hg19/chr${ii}.GL.hg38ToHg19.snp.pos |grep -v 'Unmap' |awk '{if($1==$5) print $1,$2,$3}' >EINSTEIN_hg19/chr${ii}.GL.hg38.snp.pos.bed && ## snp pos in hg38
	bedtools getfasta -fi human_ancestor_GRCh37_e59/human_ancestor_${ii}.fa -bed EINSTEIN_hg19/chr${ii}.GL.hg19.snp.pos.bed -tab |sed 's/[:,-]/\t/g' > EINSTEIN_AA/chr${ii}.GL.hg19.snp.AA.bed && ## snp pos in hg19 with aa
	paste EINSTEIN_hg19/chr${ii}.GL.hg38.snp.pos.bed EINSTEIN_AA/chr${ii}.GL.hg19.snp.AA.bed |awk -v OFS='\t' '{print $1,$3,$7}' |sed 's/[a-z]/\U&/g' |bgzip >EINSTEIN_AA/chr${ii}.GL.hg38.snp.AA.gz && ##snp pos in hg38 wuth aa
	tabix -s1 -b2 -e2 EINSTEIN_AA/chr${ii}.GL.hg38.snp.AA.gz &&
	bcftools annotate -a EINSTEIN_AA/chr${ii}.GL.hg38.snp.AA.gz -c CHROM,POS,INFO/AA -h header.AA.txt EINSTEIN_filter/chr${ii}.GL.hg38.snp.vcf.gz -Oz -o EINSTEIN_withAA/chr${ii}.GL.hg38.snp.withAA.vcf.gz && ## filtered file with aa
	rm EINSTEIN_filter/chr${ii}.GL.hg38.snp.vcf.gz */chr${ii}.GL.*bed EINSTEIN_hg19/chr${ii}.GL.hg38ToHg19.snp.pos &&
	bcftools view -S LGP_Proband.2nd_unrelated.list EINSTEIN_withAA/chr${ii}.GL.hg38.snp.withAA.vcf.gz -Oz -o EINSTEIN_withAA/Proband/LGP_Proband.${ii}.GL.hg38.snp.withAA.vcf.gz > EINSTEIN_withAA/Proband/LGP_Proband.${ii}.GL.hg38.snp.withAA.vcf.log 2>&1 && ## population 1
	bcftools view -S LGP_Control_Only.2nd_unrelated.list EINSTEIN_withAA/chr${ii}.GL.hg38.snp.withAA.vcf.gz -Oz -o EINSTEIN_withAA/Control/LGP_Control_Only.${ii}.GL.hg38.snp.withAA.vcf.gz > EINSTEIN_withAA/Control/LGP_Control_Only.${ii}.GL.hg38.snp.withAA.vcf.log 2>&1 && ## population 2
	bcftools index EINSTEIN_withAA/Proband/LGP_Proband.${ii}.GL.hg38.snp.withAA.vcf.gz &&
	bcftools index EINSTEIN_withAA/Control/LGP_Control_Only.${ii}.GL.hg38.snp.withAA.vcf.gz
}&
done
wait 

## 合并染色体
bcftools concat EINSTEIN_withAA/Proband/LGP_Proband.{1..22}.GL.hg38.snp.withAA.vcf.gz -Oz -o EINSTEIN_withAA/Proband/LGP_Proband.GL.hg38.snp.withAA.vcf.gz
bcftools concat EINSTEIN_withAA/Control/LGP_Control_Only.{1..22}.GL.hg38.snp.withAA.vcf.gz -Oz -o EINSTEIN_withAA/Control/LGP_Control.GL.hg38.snp.withAA.vcf.gz
cp LGP_Proband.2nd_unrelated.list EINSTEIN_withAA/Proband/LGP_Proband.GL.hg38.snp.withAA.list
cp LGP_Control_Only.2nd_unrelated.list EINSTEIN_withAA/Control/LGP_Control.GL.hg38.snp.withAA.list

echo "######################################################### check the EINSTEIN_withAA/Proband/LGP_Proband.GL.hg38.snp.withAA.vcf.gz ##############################################"
echo "######################################################### if ok, try: ##############################################"
echo "##### nohup bash command_list.sh EINSTEIN_withAA/Proband LGP_Proband.GL.hg38.snp.withAA output_Proband > output_Proband/LGP_Proband.GL.hg38.snp.withAA.log log 2>&1 #####"
echo "##### nohup bash command_list.sh EINSTEIN_withAA/Control LGP_Control.GL.hg38.snp.withAA output_Control > output_Control/LGP_Control.GL.hg38.snp.withAA.log log 2>&1 #####"

#nohup bash command_list.sh EINSTEIN_withAA/Proband LGP_Proband.GL.hg38.snp.withAA output_Proband > output_Proband/LGP_Proband.GL.hg38.snp.withAA.log log 2>&1 
#nohup bash command_list.sh EINSTEIN_withAA/Control LGP_Control.GL.hg38.snp.withAA output_Control > output_Control/LGP_Control.GL.hg38.snp.withAA.log log 2>&1 

