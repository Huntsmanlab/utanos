#!/bin/bash

# batch_script_qdnaseq.sh
# Description: Batch script for for creating the qdnaseq relative copy numbers

# Temporary files hierarchy:
# copy_number_calling
# .
# ├── runscripts
# │   ├── runQDNAseq.XXX.R
# │   ├── runQDNAseq.YYY.R
# │   └── ...
# ├── logs
# │   ├── runQDNAseq.XXX.R
# │   ├── runQDNAseq.YYY.R
# │   └── ...

# Permanent files hierarchy:
# copy_number_calling
# .
# ├── qdnaseq
# │   ├── XXkb_copyNumbersSegmented.rds
# │   ├── XXkb_copyNumbersSegmented.tsv
# │   ├── XXkb_copyNumbersSmooth.tsv
# │   └── ...
# ├── batch_script.sh
# └── rascal


# All paths are absolute for now
# Important main paths to set
PROJECT_DIR='copy_number_calling'
TMP_PATH=/projects/molonc/scratch/madouglas/cn_signatures/${PROJECT_DIR}
OUTPUT_PATH=/projects/molonc/huntsman_lab/madouglas/cn_signatures_shallowWGS/manual_run/${PROJECT_DIR}
BAMS_PATH=/projects/molonc/huntsman_lab/madouglas/cn_signatures_shallowWGS/manual_run/aligned

# Paths to make
RUNSCRIPTS_PATH=${TMP_PATH}/runscripts
LOGS_PATH=${TMP_PATH}/logs
mkdir -p $RUNSCRIPTS_PATH $LOGS_PATH

ls ${BAMS_PATH}/*/*.pe.bwa.sorted.dup_rm.bam > ${TMP_PATH}/bams_input.se.bams.list.txt

tester=$@
for binsize in $tester
do
  echo "binsize: $binsize"

	echo "#!/bin/bash
#SBATCH --partition=upgrade
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=${LOGS_PATH}/${binsize}.std_out.txt
#SBATCH --job-name=cn_signatures_CN_calling
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --workdir=${TMP_PATH}

Rscript_4.0 runQDNAseq.R ${TMP_PATH}/bams_input.se.bams.list.txt \
                         ${OUTPUT_PATH}/qdnaseq \
                         ${binsize}
" > ${RUNSCRIPTS_PATH}/${binsize}.qdnaseq.run.job
	echo "Submitting ${RUNSCRIPTS_PATH}/${binsize}.qdnaseq.run.job to the cluster"
	sbatch ${RUNSCRIPTS_PATH}/${binsize}.qdnaseq.run.job

done
