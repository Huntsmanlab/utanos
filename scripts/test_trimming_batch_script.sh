#!/bin/bash

# test_trimming.sh
# Description: Batch script for sWGS testing of run params
# We assume running this script from the run directory declared for SLURM
PROJECT_DIR='cn_signatures'
INPUT_PATH=/projects/molonc/huntsman_lab/madouglas/cn_signatures_shallowWGS/rawdata
OUTPUT_PATH=/projects/molonc/scratch/madouglas/${PROJECT_DIR}
PE_ADAPTERS=/projects/molonc/huntsman_lab/madouglas/bin/Trimmomatic-0.39/adapters/TruSeq3-PE.fa

mkdir -p test_trimming_logs runscripts fastqc

IFS=$'\t'
while read sow library index_i7_i5 index_i5RC flowcell lane original_source alert_code alert_notes fastq_path_1 fastq_path_2 index_i5 index_i7 i7_seq i5_seq_hiseq2500 sample_id batch status tumour_type qubit_ng_uL uL_for_10ng avg_bp_from_agilent
do
	FORFILE=$(echo 'HFFVTCCX2_3_1_'$index_i7_i5'_150bp.concat.fastq.gz')
	REVFILE=$(echo 'HFFVTCCX2_3_2_'$index_i7_i5'_150bp.concat.fastq.gz')
	mkdir -p test_trimming/$sample_id
  mkdir -p fastqc/type1 fastqc/type2 fastqc/type3

	echo "#!/bin/bash

# Change the keepBothReads param for adapters
java -jar /projects/molonc/huntsman_lab/madouglas/bin/Trimmomatic-0.39/trimmomatic-0.39.jar \
  PE -threads 16 -phred33 \
  $INPUT_PATH/$library/150bp/$FORFILE \
  $INPUT_PATH/$library/150bp/$REVFILE \
  $OUTPUT_PATH/test_trimming/$sample_id/$sample_id.f.paired.fq \
  $OUTPUT_PATH/test_trimming/$sample_id/$sample_id.f.unpaired.fq \
  $OUTPUT_PATH/test_trimming/$sample_id/$sample_id.r.paired.fq \
  $OUTPUT_PATH/test_trimming/$sample_id/$sample_id.r.unpaired.fq \
  ILLUMINACLIP:$PE_ADAPTERS:2:30:10:2:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:100

echo "${sample_id}.f.paired.fq"
cat test_trimming/$sample_id/$sample_id.f.paired.fq | awk '{if(NR%4==2) print length($1)}' > $sample_id.input.readslength.txt
textHistogram -binSize=10 $sample_id.input.readslength.txt
echo "${sample_id}.r.paired.fq"
cat test_trimming/$sample_id/$sample_id.r.paired.fq | awk '{if(NR%4==2) print length($1)}' > $sample_id.input.readslength.txt
textHistogram -binSize=10 $sample_id.input.readslength.txt
echo "${sample_id}.f.unpaired.fq"
cat test_trimming/$sample_id/$sample_id.f.unpaired.fq | awk '{if(NR%4==2) print length($1)}' > $sample_id.input.readslength.txt
textHistogram -binSize=10 $sample_id.input.readslength.txt
echo "${sample_id}.r.unpaired.fq"
cat test_trimming/$sample_id/$sample_id.r.unpaired.fq | awk '{if(NR%4==2) print length($1)}' > $sample_id.input.readslength.txt
textHistogram -binSize=10 $sample_id.input.readslength.txt

fastqc $OUTPUT_PATH/test_trimming/$sample_id/$sample_id.f.paired.fq -o fastqc/type1/


# Change the keepBothReads param for adapters & use MAXINFO rather than MINLEN
java -jar /projects/molonc/huntsman_lab/madouglas/bin/Trimmomatic-0.39/trimmomatic-0.39.jar \
	PE -threads 16 -phred33 \
	$INPUT_PATH/$library/150bp/$FORFILE \
	$INPUT_PATH/$library/150bp/$REVFILE \
	$OUTPUT_PATH/test_trimming/$sample_id/$sample_id.f.paired.fq \
	$OUTPUT_PATH/test_trimming/$sample_id/$sample_id.f.unpaired.fq \
	$OUTPUT_PATH/test_trimming/$sample_id/$sample_id.r.paired.fq \
	$OUTPUT_PATH/test_trimming/$sample_id/$sample_id.r.unpaired.fq \
	ILLUMINACLIP:$PE_ADAPTERS:2:30:10:2:false SLIDINGWINDOW:4:15 MAXINFO:100:0.2

echo "${sample_id}.f.paired.fq"
cat test_trimming/$sample_id/$sample_id.f.paired.fq | awk '{if(NR%4==2) print length($1)}' > $sample_id.input.readslength.txt
textHistogram -binSize=10 $sample_id.input.readslength.txt
echo "${sample_id}.r.paired.fq"
cat test_trimming/$sample_id/$sample_id.r.paired.fq | awk '{if(NR%4==2) print length($1)}' > $sample_id.input.readslength.txt
textHistogram -binSize=10 $sample_id.input.readslength.txt
echo "${sample_id}.f.unpaired.fq"
cat test_trimming/$sample_id/$sample_id.f.unpaired.fq | awk '{if(NR%4==2) print length($1)}' > $sample_id.input.readslength.txt
textHistogram -binSize=10 $sample_id.input.readslength.txt
echo "${sample_id}.r.unpaired.fq"
cat test_trimming/$sample_id/$sample_id.r.unpaired.fq | awk '{if(NR%4==2) print length($1)}' > $sample_id.input.readslength.txt
textHistogram -binSize=10 $sample_id.input.readslength.txt

fastqc $OUTPUT_PATH/test_trimming/$sample_id/$sample_id.f.paired.fq -o fastqc/type2/

# Change the keepBothReads param for adapters & use MAXINFO rather than MINLEN
java -jar /projects/molonc/huntsman_lab/madouglas/bin/Trimmomatic-0.39/trimmomatic-0.39.jar \
	PE -threads 16 -phred33 \
	$INPUT_PATH/$library/150bp/$FORFILE \
	$INPUT_PATH/$library/150bp/$REVFILE \
	$OUTPUT_PATH/test_trimming/$sample_id/$sample_id.f.paired.fq \
	$OUTPUT_PATH/test_trimming/$sample_id/$sample_id.f.unpaired.fq \
	$OUTPUT_PATH/test_trimming/$sample_id/$sample_id.r.paired.fq \
	$OUTPUT_PATH/test_trimming/$sample_id/$sample_id.r.unpaired.fq \
	ILLUMINACLIP:$PE_ADAPTERS:2:30:10:2:false SLIDINGWINDOW:4:15 MAXINFO:100:0.5

echo "${sample_id}.f.paired.fq"
cat test_trimming/$sample_id/$sample_id.f.paired.fq | awk '{if(NR%4==2) print length($1)}' > $sample_id.input.readslength.txt
textHistogram -binSize=10 $sample_id.input.readslength.txt
echo "${sample_id}.r.paired.fq"
cat test_trimming/$sample_id/$sample_id.r.paired.fq | awk '{if(NR%4==2) print length($1)}' > $sample_id.input.readslength.txt
textHistogram -binSize=10 $sample_id.input.readslength.txt
echo "${sample_id}.f.unpaired.fq"
cat test_trimming/$sample_id/$sample_id.f.unpaired.fq | awk '{if(NR%4==2) print length($1)}' > $sample_id.input.readslength.txt
textHistogram -binSize=10 $sample_id.input.readslength.txt
echo "${sample_id}.r.unpaired.fq"
cat test_trimming/$sample_id/$sample_id.r.unpaired.fq | awk '{if(NR%4==2) print length($1)}' > $sample_id.input.readslength.txt
textHistogram -binSize=10 $sample_id.input.readslength.txt

fastqc $OUTPUT_PATH/test_trimming/$sample_id/$sample_id.f.paired.fq -o fastqc/type3/

  " > runscripts/$sample_id.analysis.bash

  	echo "#!/bin/bash
  #SBATCH --partition=upgrade
  #SBATCH --time=02:00:00
  #SBATCH --nodes=1
  #SBATCH --ntasks=1
  #SBATCH --output=${OUTPUT_PATH}/test_trimming_logs/$sample_id.std_out_1.txt
  #SBATCH --job-name=trimming_testing
  #SBATCH --mem=32G
  #SBATCH --cpus-per-task=20
  #SBATCH --workdir=${OUTPUT_PATH}
  sh runscripts/$sample_id.analysis.bash" > runscripts/$sample_id.run.job
  		echo "Submitting runscripts/$sample_id.run.job to the cluster"
  		sbatch runscripts/$sample_id.run.job

  done < "metadata_no_colnames_tester.tsv"
