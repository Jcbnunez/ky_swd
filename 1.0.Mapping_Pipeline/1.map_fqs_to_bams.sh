#!/usr/bin/env bash
#
#SBATCH -J map # A single job name for the array
#SBATCH -c 1
#SBATCH -N 1 # on one node
#SBATCH -t 12:00:00 #<= this may depend on your resources
#SBATCH --mem 80G #<= this may depend on your resources
#SBATCH -p bluemoon
#SBATCH -o ./slurmOutput/map.%A_%a.out # Standard output
#SBATCH --array=1

###-20


## sanity check post load
my_job_header
####
####
####
bwamem2=/netfiles/nunezlab/Shared_Resources/Software/bwa-mem2-2.2.1_x64-linux/bwa-mem2.avx2
samtools=/netfiles/nunezlab/Shared_Resources/Software/samtools-1.19/samtools

#### user inputs
metadata="/gpfs2/scratch/jcnunez/Dsu.prelim.data/SampleId_metadata.txt"
seq_id=$( cat $metadata  | sed '1d' | awk '{print $2}' | sed "${SLURM_ARRAY_TASK_ID}q;d" )
QUERY_NAME=$seq_id

mkdir -p bam_files

PATH_TO_ASSEMBLY="/netfiles/nunezlab/D_suzukii_resources/genomes/MG_INRA2024_genome/dsu_isojap1_chrlevel_ncbi.fa"
dir_data="/gpfs2/scratch/jcnunez/Dsu.prelim.data/Clean_Reads"

NTHREADS=$SLURM_CPUS_ON_NODE
SAMTOOLS_THREADS=$(($NTHREADS-1))

RGLB='LIB-'$QUERY_NAME
RGSM=$QUERY_NAME

#####
#####
#####

#mapping, on filtre q20, puis on sorte par query_name (important pour le fixmate lui-meme necessaire pour le markdup)
$bwamem2 mem \
-M -t $NTHREADS \
-R @RG\\tID:${QUERY_NAME}\\tPL:ILLUMINA\\tLB:${RGLB}\\tSM:${RGSM} $PATH_TO_ASSEMBLY \
$dir_data/$QUERY_NAME.clean.F.gz \
$dir_data/$QUERY_NAME.clean.R.gz | \
$samtools sort -T $QUERY_NAME -n --threads $SAMTOOLS_THREADS \
-O BAM -o bam_files/$QUERY_NAME.sorted.bam -

######
mkdir -p bam_files_fm

$samtools fixmate -rpcm --threads $SAMTOOLS_THREADS \
bam_files/$QUERY_NAME.sorted.bam bam_files_fm/$QUERY_NAME.sorted.fm.bam

##clean
rm bam_files/$QUERY_NAME.sorted.bam

######
mkdir -p bam_files_sort2

$samtools sort -T $QUERY_NAME --threads $SAMTOOLS_THREADS \
bam_files_fm/$QUERY_NAME.sorted.fm.bam > bam_files_sort2/$QUERY_NAME.sorted.bam

##clean
rm bam_files_fm/$QUERY_NAME.sorted.fm.bam

######
######
mkdir -p bam_files_final


$samtools markdup --threads $SAMTOOLS_THREADS \
-r bam_files_sort2/$QUERY_NAME.sorted.bam bam_files_final/$QUERY_NAME.dedup.bam

### index
$samtools index -@ $SAMTOOLS_THREADS \
bam_files_final/$QUERY_NAME.dedup.bam

#### stats
mkdir -p bam_stats

$samtools stats --threads $SAMTOOLS_THREADS \
bam_files_final/$QUERY_NAME.dedup.bam > bam_stats/$QUERY_NAME.dedup.bamstats

rm bam_files_sort2/$QUERY_NAME.sorted.bam


