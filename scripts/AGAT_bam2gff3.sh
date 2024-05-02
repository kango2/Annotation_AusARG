#!/bin/bash
#PBS -N AGAT
#PBS -l ncpus=16,walltime=4:00:00,storage=gdata/if89+gdata/xl04,mem=100GB
#PBS -j oe


set -e
module use /g/data/if89/apps/modulefiles
module unload AGAT samtools perllib parallel
module load AGAT/1.4.0 parallel/20191022


# bam to gff3 conversion using AGAT
ls ${input}/*.bam | parallel --jobs ${PBS_NCPUS} '
  export base=$(basename {} .bam)
  agat_convert_minimap2_bam2gff.pl -i {} -o ${output}/${base}.gff3
'
