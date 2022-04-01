#!/bin/bash

# batch_script.sh
# Description: Batch script for sWGS data processing

# All paths are absolute for now
# Important main paths to set
PROJECT_DIR='cn_signatures_shallowWGS'
TMP_PATH=/projects/molonc/scratch/madouglas/cn_signatures
OUTPUT_PATH=/projects/molonc/huntsman_lab/madouglas/${PROJECT_DIR}
REF_GENOME=${OUTPUT_PATH}/../bin/reference_genomes/dlp_refdata/human/GRCh37-lite.fa
SE_ADAPTERS=${OUTPUT_PATH}/../bin/Trimmomatic-0.39/adapters/TruSeq3-SE.fa
PE_ADAPTERS=${OUTPUT_PATH}/../bin/Trimmomatic-0.39/adapters/TruSeq3-PE.fa


IFS=$'\t'
while read sample_id batch repeat_run tissue sample_type status hist cancer_type grade sequence library fastq_path_1 fastq_path_2
do
	# if [[ $batch == 1 ]] || [[ $batch == 2 ]] || [[ $batch == 3 ]] || \
	# 	[[ $batch == 4 ]] || [[ $batch == 5 ]] || [[ $batch == 6 ]] || \
	# 	[[ $batch == 7 ]] || [[ $batch == 8 ]] || [[ $batch == 9 ]] || \
	# 	[[ $batch == 10 ]] || [[ $batch == 11 ]];
	if [[ $repeat_run == 1 ]];
	then

		# Make pipeline directories
		RUNSCRIPTS=${TMP_PATH}/runscripts/${sample_id}
		TRIMMED=${TMP_PATH}/trimmed/${sample_id}
		ALIGNED=${TMP_PATH}/aligned/${sample_id}

		FASTQC_RAW=${OUTPUT_PATH}/manual_run/fastQC/${sample_id}/raw
		FASTQC_PROC=${OUTPUT_PATH}/manual_run/fastQC/${sample_id}/post_alignment
		LOGS=${TMP_PATH}/logs/${sample_id}

		mkdir -p $RUNSCRIPTS $TRIMMED $ALIGNED
		mkdir -p $LOGS $FASTQC_RAW $FASTQC_PROC

		# echo $fastq_path_1
		echo "#!/bin/bash
source /projects/molonc/huntsman_lab/madouglas/bin/.conda/envs/rnaseq_env/bin/activate

fastqc $fastq_path_1 -o $FASTQC_RAW/
fastqc $fastq_path_2 -o $FASTQC_RAW/

java -jar /projects/molonc/huntsman_lab/madouglas/bin/Trimmomatic-0.39/trimmomatic-0.39.jar \
	PE -threads 16 -phred33 \
	$fastq_path_1 \
	$fastq_path_2 \
	$TRIMMED/$sample_id.f.paired.fq $TRIMMED/$sample_id.f.unpaired.fq \
	$TRIMMED/$sample_id.r.paired.fq $TRIMMED/$sample_id.r.unpaired.fq \
	ILLUMINACLIP:$PE_ADAPTERS:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40 MAXINFO:100:0.5

echo "${sample_id}.f.paired.fq"
cat $TRIMMED/$sample_id.f.paired.fq | awk '{if(NR%4==2) print length($1)}' > $TRIMMED/$sample_id.input.readslength.txt
textHistogram -binSize=10 $TRIMMED/$sample_id.input.readslength.txt
echo "${sample_id}.f.unpaired.fq"
cat $TRIMMED/$sample_id.f.unpaired.fq | awk '{if(NR%4==2) print length($1)}' > $TRIMMED/$sample_id.input.readslength.txt
textHistogram -binSize=10 $TRIMMED/$sample_id.input.readslength.txt
echo "${sample_id}.r.paired.fq"
cat $TRIMMED/$sample_id.r.paired.fq | awk '{if(NR%4==2) print length($1)}' > $TRIMMED/$sample_id.input.readslength.txt
textHistogram -binSize=10 $TRIMMED/$sample_id.input.readslength.txt
echo "${sample_id}.r.unpaired.fq"
cat $TRIMMED/$sample_id.r.unpaired.fq | awk '{if(NR%4==2) print length($1)}' > $TRIMMED/$sample_id.input.readslength.txt
textHistogram -binSize=10 $TRIMMED/$sample_id.input.readslength.txt

#### Paired-End Alignment
bwa-mem2 mem -M -t 24 $REF_GENOME \
		$TRIMMED/$sample_id.f.paired.fq \
		$TRIMMED/$sample_id.r.paired.fq \
		> $ALIGNED/$sample_id.pe.bwa.sam
samtools view -h -b -S $ALIGNED/$sample_id.pe.bwa.sam > $ALIGNED/$sample_id.pe.bwa.bam
samtools sort -m 2G -@ 24 -o $ALIGNED/$sample_id.pe.bwa.sorted.bam $ALIGNED/$sample_id.pe.bwa.bam
samtools index -@ 24 $ALIGNED/$sample_id.pe.bwa.sorted.bam
java -jar /projects/molonc/huntsman_lab/madouglas/bin/picard.jar MarkDuplicates \
		I=$ALIGNED/$sample_id.pe.bwa.sorted.bam \
		O=$ALIGNED/$sample_id.pe.bwa.sorted.dup_rm.bam \
		M=$ALIGNED/$sample_id.pe.marked_dup_metrics.txt
samtools index -@ 24 $ALIGNED/$sample_id.pe.bwa.sorted.dup_rm.bam
fastqc $ALIGNED/$sample_id.pe.bwa.sorted.dup_rm.bam -o $FASTQC_PROC/
samtools coverage $ALIGNED/$sample_id.pe.bwa.sorted.dup_rm.bam > $OUTPUT_PATH/manual_run/coverage/$sample_id.pe.coverageTable.tsv

#### Single-End Alignment (forward strand)
bwa-mem2 mem -M -t 24 $REF_GENOME \
		$TRIMMED/$sample_id.f.unpaired.fq \
		> $ALIGNED/$sample_id.se.bwa.sam
samtools view -h -b -S $ALIGNED/$sample_id.se.bwa.sam > $ALIGNED/$sample_id.se.bwa.bam
samtools sort -m 2G -@ 24 -o $ALIGNED/$sample_id.se.bwa.sorted.bam $ALIGNED/$sample_id.se.bwa.bam
samtools index -@ 24 $ALIGNED/$sample_id.se.bwa.sorted.bam
# Grab forward strand from P.E. post_alignment
samtools view -F 16 -o $ALIGNED/$sample_id.pe.f.bwa.sorted.bam $ALIGNED/$sample_id.pe.bwa.sorted.bam
samtools merge -f -@ 24 $ALIGNED/$sample_id.se.merged.bam \
		$ALIGNED/$sample_id.se.bwa.sorted.bam \
		$ALIGNED/$sample_id.pe.f.bwa.sorted.bam
# Re-Sort now merged BAM
samtools sort -m 2G -@ 24 -o $ALIGNED/$sample_id.se.merged.sorted.bam $ALIGNED/$sample_id.se.merged.bam
samtools index -@ 24 $ALIGNED/$sample_id.se.merged.sorted.bam
java -jar /projects/molonc/huntsman_lab/madouglas/bin/picard.jar MarkDuplicates \
		I=$ALIGNED/$sample_id.se.merged.sorted.bam \
		O=$ALIGNED/$sample_id.se.merged.sorted.dup_rm.bam \
		M=$ALIGNED/$sample_id.se.marked_dup_metrics.txt
samtools index -@ 24 $ALIGNED/$sample_id.se.merged.sorted.dup_rm.bam
fastqc $ALIGNED/$sample_id.se.merged.sorted.dup_rm.bam -o $FASTQC_PROC/
samtools coverage $ALIGNED/$sample_id.se.merged.sorted.dup_rm.bam > $OUTPUT_PATH/manual_run/coverage/$sample_id.se.coverageTable.tsv

#### Clean-up
rm $ALIGNED/$sample_id.pe.bwa.sam
rm $ALIGNED/$sample_id.pe.bwa.bam
rm $ALIGNED/$sample_id.pe.bwa.sorted.bam
rm $ALIGNED/$sample_id.pe.bwa.sorted.bam.bai

rm $ALIGNED/$sample_id.se.bwa.sam
rm $ALIGNED/$sample_id.se.bwa.bam
rm $ALIGNED/$sample_id.se.bwa.sorted.bam
rm $ALIGNED/$sample_id.se.bwa.sorted.bam.bai
rm $ALIGNED/$sample_id.pe.f.bwa.sorted.bam
rm $ALIGNED/$sample_id.se.merged.bam
rm $ALIGNED/$sample_id.se.merged.sorted.bam
rm $ALIGNED/$sample_id.se.merged.sorted.bam.bai

" > $RUNSCRIPTS/${sample_id}.analysis.bash

		echo "#!/bin/bash
#SBATCH --partition=upgrade
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=${LOGS}/${sample_id}.std_out.txt
#SBATCH --job-name=cn_signatures_alignment
#SBATCH --mem=128G
#SBATCH --cpus-per-task=32
#SBATCH --workdir=${TMP_PATH}
sh runscripts/$sample_id/${sample_id}.analysis.bash" > $RUNSCRIPTS/${sample_id}.run.job
		echo "Submitting $RUNSCRIPTS/${sample_id}.run.job to the cluster"
		sbatch $RUNSCRIPTS/${sample_id}.run.job
	fi
done < "/projects/molonc/huntsman_lab/madouglas/cn_signatures_shallowWGS/metadata/metadata_unix.tsv"
