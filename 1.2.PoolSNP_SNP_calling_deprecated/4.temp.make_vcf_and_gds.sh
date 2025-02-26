#!/usr/bin/env bash
#
#SBATCH -J dockerMap # A single job name for the array
#SBATCH -c 1
#SBATCH -N 1 # on one node
#SBATCH -t 5:00:00 #<= this may depend on your resources
#SBATCH --mem 45G #<= this may depend on your resources
#SBATCH -o ./slurmOutput/RunDest.%A_%a.out # Standard output
#SBATCH -p bluemoon
#!/usr/bin/env bash


#module load htslib bcftools intel/18.0 intelmpi/18.0 parallel R/3.6.3
module load Rtidyverse

### Load picard
#SNPEFF=/users/j/n/jnunez2/software/snpEff/exec/snpeff
PICARD=/netfiles/nunezlab/Shared_Resources/Software/picard/build/libs/picard-2.27.5-5-gbdfea14-SNAPSHOT-all.jar
JAVAMEM=49G
PROJNAME=SWD_KY
DICT=/netfiles/nunezlab/D_suzukii_resources/genomes/MG_INRA2024_genome/dsu_isojap1_chrlevel_ncbi.dict
vcf2gds=/gpfs2/scratch/jcnunez/test_DEST/poolSNP/snpCalling/innerscripts/vcf2gds.R
tabix=/netfiles/nunezlab/Shared_Resources/Software/htslib/tabix
bgzip=/netfiles/nunezlab/Shared_Resources/Software/htslib/bgzip

###
##### params
  popSet="PoolSeq"
  method="PoolSNP"
  maf="001"
  mac=5
  version=ORCC
  script_dir="/gpfs2/scratch/jcnunez/test_DEST/poolSNP/snpCalling/innerscripts"
  wd="/users/j/c/jcnunez/scratch/Dsu.prelim.data/scatter"

syncPath1orig="/gpfs2/scratch/jcnunez/Dsu.prelim.data/BAM_to_SYNC/synfiles_KY/*/*masked.sync.gz"  
syncPath2orig="/gpfs2/scratch/jcnunez/Dsu.prelim.data/BAM_to_SYNC/synfiles_KY/*/*masked.sync.gz"

#####

java -Xmx$JAVAMEM \
 -jar $PICARD MergeVcfs \
  I=${wd}/sub_bcf/dest.chr2L.${popSet}.${method}.${maf}.${mac}.${version}.vcf \
  I=${wd}/sub_bcf/dest.chr2R.${popSet}.${method}.${maf}.${mac}.${version}.vcf \
  I=${wd}/sub_bcf/dest.chr3.${popSet}.${method}.${maf}.${mac}.${version}.vcf \
  I=${wd}/sub_bcf/dest.chr4.${popSet}.${method}.${maf}.${mac}.${version}.vcf \
  I=${wd}/sub_bcf/dest.chrX.${popSet}.${method}.${maf}.${mac}.${version}.vcf \
  I=${wd}/sub_bcf/dest.chrY.${popSet}.${method}.${maf}.${mac}.${version}.vcf \
  O=${wd}/$PROJNAME.${popSet}.${method}.${maf}.${mac}.${version}.vcf \
  D=$DICT

echo "fix header" #this is now fixed in PoolSNP.py
  sed -i '0,/CHROM/{s/AF,Number=1/AF,Number=A/}' ${wd}/$PROJNAME.${popSet}.${method}.${maf}.${mac}.${version}.vcf
  sed -i '0,/CHROM/{s/AC,Number=1/AC,Number=A/}' ${wd}/$PROJNAME.${popSet}.${method}.${maf}.${mac}.${version}.vcf
  sed -i '0,/CHROM/{s/AD,Number=1/AD,Number=A/}' ${wd}/$PROJNAME.${popSet}.${method}.${maf}.${mac}.${version}.vcf
  sed -i '0,/CHROM/{s/FREQ,Number=1/FREQ,Number=A/}' ${wd}/$PROJNAME.${popSet}.${method}.${maf}.${mac}.${version}.vcf

####
echo "make GDS"
  Rscript --vanilla \
  $vcf2gds \
  ${wd}/$PROJNAME.${popSet}.${method}.${maf}.${mac}.${version}.vcf

echo "bgzip & tabix"
  $bgzip -c  ${wd}/$PROJNAME.${popSet}.${method}.${maf}.${mac}.${version}.vcf > ${wd}/$PROJNAME.${popSet}.${method}.${maf}.${mac}.${version}.vcf.gz
  $tabix -p vcf ${wd}/$PROJNAME.${popSet}.${method}.${maf}.${mac}.${version}.vcf.gz
