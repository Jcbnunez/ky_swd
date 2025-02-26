#!/usr/bin/env bash
#
#SBATCH -J dockerMap # A single job name for the array
#SBATCH -c 1
#SBATCH -N 1 # on one node
#SBATCH -t 5:00:00 #<= this may depend on your resources
#SBATCH --mem 45G #<= this may depend on your resources
#SBATCH -o ./slurmOutput/RunDest.%A_%a.out # Standard output
#SBATCH -e ./slurmOutput/RunDest.%A_%a.err # Standard error
#SBATCH -p bluemoon
#SBATCH --array=1-6

#module load spack/spack-0.18.1
#spack load bcftools@1.14

PICARD=/netfiles/nunezlab/Shared_Resources/Software/picard/build/libs/picard-2.27.5-5-gbdfea14-SNAPSHOT-all.jar

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
chr_file=/gpfs2/scratch/jcnunez/Dsu.prelim.data/scatter/chr_file.txt
DICT=/netfiles/nunezlab/D_suzukii_resources/genomes/MG_INRA2024_genome/dsu_isojap1_chrlevel_ncbi.dict


##### SCRIPT ITERATES HERE ####
chr=$( cat $chr_file | sed "${SLURM_ARRAY_TASK_ID}q;d" )
echo $SLURM_ARRAY_TASK_ID
echo $chr
#####

echo "Chromosome: $chr"
date

bcf_outdir="${wd}/sub_bcf"
if [ ! -d $bcf_outdir ]; then
    mkdir $bcf_outdir
fi

outdir=$wd/sub_vcfs
cd ${wd}

ls -d ${outdir}/*.${popSet}.${method}.${maf}.${mac}.${version}.vcf.gz | sort -t"_" -k2,2 -k4g,4  | \
grep /${chr}_ > $outdir/vcfs_order.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.sort.list

echo "Concatenating"

### Requires dictionary
java -jar $PICARD MergeVcfs \
          I=$outdir/vcfs_order.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.sort.list  \
          O=$bcf_outdir/dest.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.vcf \
          D=$DICT
          
echo "done"
date