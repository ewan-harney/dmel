#!/bin/sh

#SBATCH --partition=haswell
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --mail-type=START,END,FAIL
#SBATCH -o logs/exprmnt_%a.out # Standard output
#SBATCH -e logs/exprmnt_%a.err # Standard error
#SBATCH --mail-user=ewan.harney@gmail.com

module load Nextflow/22.10.5
module load Singularity/3.8.7-GCCcore-11.2.0

nextflow run nf-core/atacseq \
--input /homes/users/eharney/scratch/Nextflow_exprmnt/samplesheet_edh_exprmnt.csv \
--read_length 50 \
--outdir /homes/users/eharney/scratch/Nextflow_exprmnt \
--fasta /homes/users/eharney/scratch/Nextflow_exprmnt/Dmel_FB_r6_v46_chr_andMt.fasta \
--gtf /homes/users/eharney/scratch/Nextflow_exprmnt/Dmel_FB_r6_v46_chr_andMt.gtf \
--mito_name chrMt \
-profile singularity \
-c /homes/users/eharney/scratch/Nextflow_exprmnt/marvin.config
