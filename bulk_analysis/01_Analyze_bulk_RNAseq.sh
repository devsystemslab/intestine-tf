#!/bin/bash
#SBATCH --job-name=Analyze_bulk_RNAseq         # Job name
#SBATCH --partition=batch_cpu                 # Partition name
#SBATCH --qos=3h
#SBATCH --mem=50G                       # Total memory
#SBATCH --cpus-per-task=1                # 1 CPU
#SBATCH --output=/home/adaml9/logs/Analyze_bulk_RNAseq.out
#SBATCH --error=/home/adaml9/logs/Analyze_bulk_RNAseq.err

export PATH=/projects/site/pred/ihb-global/ihb-genomics-framework/conda/nf-core/bin:$PATH
export JAVA_CMD=/projects/site/pred/ihb-global/ihb-genomics-framework/conda/nf-core/bin/java

DATA_DIR="/projects/site/pred/ihb-intestine-evo/lukas_area/eec_fate_project/bulkrna/processed"
RESULTS_DIR="/projects/site/pred/ihb-g-deco/USERS/adaml9/tmp/bulk_results" 
WORK_DIR="/projects/site/pred/ihb-g-deco/USERS/adaml9/tmp/bulk_work"
WORKFLOW_DIR="/projects/site/pred/ihb-global/ihb-genomics-framework/workflows/nf-deseq2"

cd $DATA_DIR

# Load nextflow
source ~/.bashrc
conda activate /projects/site/pred/ihb-global/ihb-genomics-framework/conda/nf-core

nextflow run $WORKFLOW_DIR \
     --input $DATA_DIR/samplesheet.csv \
     --contrasts $DATA_DIR/contrasts.csv \
     --matrix $DATA_DIR/counts.txt \
     --selected_genes $DATA_DIR/selected_genes.txt \
     --outdir deseq2 \
     -profile singularity