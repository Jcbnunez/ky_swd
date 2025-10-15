#!/usr/bin/env bash
#
#SBATCH -J CLARKp # A single job name for the array
#SBATCH -c 1
#SBATCH -N 1 # on one node
#SBATCH -t 12:00:00 #<= this may depend on your resources
#SBATCH --mem 80G #<= this may depend on your resources
#SBATCH -p bluemoon
#SBATCH -o ./slurmOutput/CLARKp.%A_%a.out # Standard output
#SBATCH --array=2-20

my_job_header

#### Clark clean pipeline
fastp=/netfiles/nunezlab/Shared_Resources/Software/fastp_latest/fastp

clark=/netfiles/nunezlab/Shared_Resources/Software/CLARKV1.3.0.0/exe/CLARK
clarkl=/netfiles/nunezlab/Shared_Resources/Software/CLARKV1.3.0.0/exe/CLARK-l
#clarks=/netfiles/nunezlab/Shared_Resources/Software/CLARKV1.3.0.0/exe/CLARK-S

####
#summary_csv.awk
parseclarkcsv=/netfiles/nunezlab/Shared_Resources/Software/CLARKV1.3.0.0/summary_csv.awk
####

#### User given data
## sanity check post load
####

#### user inputs
metadata="/gpfs2/scratch/jcnunez/Dsu.prelim.data/SampleId_metadata.txt"
seq_id=$( cat $metadata  | sed '1d' | awk '{print $2}' | sed "${SLURM_ARRAY_TASK_ID}q;d" )
echo $seq_id
dirfq=/gpfs2/scratch/jcnunez/Dsu.prelim.data/Reads_Merged #assumed to be coded as ${id}.fq_1.gz and ${id}.fq_2.gz for read1 and read2 data
TARGET=/gpfs2/scratch/jcnunez/Dsu.prelim.data/droso_clark/target_clean.txt
DB_folder=/gpfs2/scratch/jcnunez/Dsu.prelim.data/droso_clark/droso_clean_db
id=$seq_id #ID prefix to provide when launching the shel script

####
#### data handling

mkdir -p res_clarkl #create a directory named res_clarkl (if it doesn't exist) to store calrkl result
mkdir -p res_clark  #create a directory named res_clark (if it doesn't exist) to store calrk result

mkdir -p res_clarkl_jsons 
mkdir -p res_clark_jsons  

mkdir -p processed_fastas 
mkdir -p res_clark/summaries
mkdir -p res_clarkl/summaries

date
#####


####################
#cleanning seqeunce: stdout option to obtain interleaved format further transformed into fasta using awk one liner
####################

$fastp -i $dirfq"/"$id'.F.merged.fastq.gz' -I $dirfq"/"$id'.R.merged.fastq.gz' \
--stdout --merge --include_unmerged -h $id'.html' \
-j $id'.json' | awk '{if(NR%4==2){nn++;{print ">s"nn"\n"$0}}}' - > processed_fastas"/"$id'.fasta'

date


######################
##CLARK-l analysis
####################


$clarkl -T $TARGET -D $DB_folder -O processed_fastas"/"$id'.fasta' -R $id -n 1 -m 0 -s 2

date

mv ${id}.csv res_clarkl/
mv ${id}.html res_clarkl_jsons
mv ${id}.json res_clarkl_jsons

######################
##CLARK analysis
####################

$clark -T $TARGET -D $DB_folder -O processed_fastas"/"$id'.fasta' -R $id -n 1 -m 0 -s 2

date

mv ${id}.csv res_clark/
mv ${id}.html res_clark_jsons
mv ${id}.json res_clark_jsons

#### Clean house

rm processed_fastas"/"$id'.fasta'

###########
###summarize statistics (using awk script)
############

nkmin='1 5'
conf='0.9 0.95'
resdir='res_clark/ res_clarkl/'

for dd in $resdir
do
cd ./${dd}
 for ii in $nkmin
  do
  for jj in $conf
   do
    gawk -f $parseclarkcsv -v nmin_kmer=$ii -v conf_thr=$jj $id.csv > ./summaries/$id"_nkmin"$ii"_conf"$jj".summary"
   done
 done
gzip ${id}.csv
cd ../
done

date
