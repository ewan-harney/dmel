#!/bin/sh

STRAIN=${1}
READ_LENGTH=${2}

module load RepeatMasker/4.0.7-foss-2016b
module load SHRiMP/2.2.3
module load BLAT/3.5-foss-2016b
module load phrap/1.09
module load fastagrep/2.0
module load BWA/0.7.15
module load SAMtools/1.4-foss-2016b
module load BCFtools/1.3.1
module load Perl/5.22.1-foss-2016b
module load FASTX-Toolkit/0.0.14-foss-2016b

echo 'Read length is ' $READ_LENGTH

perl /gpfs42/projects/lab_jgonzalez/gonzalez_lab/programs/Tlex/tlex-open-v3.0.pl \
-O $STRAIN \
-A $READ_LENGTH \
-pairends yes \
-s drosophila \
-noFilterTE \
-T /gpfs42/robbyfs/scratch/lab_jgonzalez/eharney/Tlex/input/TElist_2417TEs_ISO1_nochr.txt \
-M /gpfs42/robbyfs/scratch/lab_jgonzalez/eharney/Tlex/input/TEcopies_2417TEs_ISO1_nochr.txt \
-G /gpfs42/robbyfs/scratch/lab_jgonzalez/eharney/genomes/Flybase_Dmel_r6-46/Dmel_FB_r6_v46_nochr.fasta \
-R /gpfs42/robbyfs/scratch/lab_jgonzalez/eharney/Tlex/$STRAIN
