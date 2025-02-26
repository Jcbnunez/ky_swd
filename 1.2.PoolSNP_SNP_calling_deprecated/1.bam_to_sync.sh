#!/usr/bin/env bash
#
#SBATCH -J make_sync # A single job name for the array
#SBATCH -c 5
#SBATCH -N 1 # on one node
#SBATCH -t 20:00:00 #<= this may depend on your resources
#SBATCH --mem 60G #<= this may depend on your resources
#SBATCH -o ./slurmOutput/make_sync.%A_%a.out # Standard output
#SBATCH -p bluemoon
#SBATCH --array=1

##-20


##user parameters
JAVAMEM=59G
threads=5
ith=${SLURM_ARRAY_TASK_ID}
####SAMP=TEST
metadata="/gpfs2/scratch/jcnunez/Dsu.prelim.data/SampleId_metadata.txt"
seq_id=$( cat $metadata  | sed '1d' | awk '{print $2}' | sed "${SLURM_ARRAY_TASK_ID}q;d" )

## file needs
root=/gpfs2/scratch/jcnunez/Dsu.prelim.data/bam_files_final
####root=./
genome=/netfiles/nunezlab/D_suzukii_resources/genomes/MG_INRA2024_genome/dsu_isojap1_chrlevel_ncbi.fa
genome_pickled=/netfiles/nunezlab/D_suzukii_resources/genomes/MG_INRA2024_genome/dsu_isojap1_chrlevel_ncbi.fa.pickled.ref
tegff=/netfiles/nunezlab/D_suzukii_resources/repeat_masking/for_MG_INRA2024_genome/dsu_isojap1_chrlevel_ncbi.fa.out.gff


### program dependencies
gatk3=/netfiles/nunezlab/Shared_Resources/Software/gatk3/gatk/GenomeAnalysisTK.jar
samtools=/netfiles/nunezlab/Shared_Resources/Software/samtools-1.19/samtools
tabix=/netfiles/nunezlab/Shared_Resources/Software/htslib/tabix
bgzip=/netfiles/nunezlab/Shared_Resources/Software/htslib/bgzip

##need python3
module load python3.10-anaconda/2023.03-1


Mpileup2Sync=/netfiles/nunezlab/Shared_Resources/Software/DESTv2/mappingPipeline/scripts/Mpileup2Sync.py
MaskSYNC_snape_complete=/netfiles/nunezlab/Shared_Resources/Software/DESTv2/mappingPipeline/scripts/MaskSYNC_snape_complete.py

#### parameters
base_quality_threshold=25
illumina_quality_coding=1.8
minIndel=5
maxsnape=0.9
max_cov=0.95 
min_cov=10 

### <*> ####
### Begin ##
### <*> ####

#### File locations
output=synfiles_KY
sample=$seq_id
echo $sample

mkdir $output
mkdir $output/$sample


###indel creator RealignerTargetCreator...

 java -jar $gatk3 -T RealignerTargetCreator \
  -nt $threads \
  -R $genome \
  -I $root/${sample}.dedup.bam \
  -o $output/$sample/${sample}.hologenome.intervals

#### inder realigner....

  java -jar $gatk3 \
  -T IndelRealigner \
  -R $genome \
  -I $root/${sample}.dedup.bam \
  -targetIntervals $output/$sample/${sample}.hologenome.intervals \
  -o $output/$sample/${sample}.contaminated_realigned.bam
  ###rm $output/$sample/${sample}.dedup.bam*


#### run mpileup step 

  $samtools mpileup \
  $output/$sample/${sample}.contaminated_realigned.bam \
  -B \
  -Q ${base_quality_threshold} \
  -f $genome > $output/$sample/${sample}_mpileup.txt

#### Transform Mpileup to Sync

  python3 $Mpileup2Sync \
  --mpileup $output/$sample/${sample}_mpileup.txt \
  --ref $genome_pickled \
  --output $output/$sample/${sample} \
  --base-quality-threshold $base_quality_threshold \
  --coding $illumina_quality_coding \
  --minIndel $minIndel

#### prepare PoolSNP output

  python3 $MaskSYNC_snape_complete \
  --sync $output/$sample/${sample}.sync.gz \
  --output $output/$sample/${sample} \
  --indel $output/$sample/${sample}.indel \
  --coverage $output/$sample/${sample}.cov \
  --mincov $min_cov \
  --maxcov $max_cov \
  --te $tegff \
  --maxsnape $maxsnape


  mv $output/$sample/${sample}_masked.sync.gz $output/$sample/${sample}.masked.sync.gz
  gunzip $output/$sample/${sample}.masked.sync.gz
  $bgzip $output/$sample/${sample}.masked.sync
  $tabix -s 1 -b 2 -e 2 $output/$sample/${sample}.masked.sync.gz

date
echo "done"