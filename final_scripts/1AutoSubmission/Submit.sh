#!/bin/bash
export TRANSLATION_OUTPUT="/g/data/xl04/hrp561/basdurnaseq"
export workingdir="/path/to/workingdir"
export genome="/g/data/xl04/jc4878/Bassiana_publication/YAHS/genome/BASDU_HifiASM_YAHS_SUP_CurV1.1.fasta.masked"
export scriptdir="/g/data/xl04/jc4878/Annotation_AusARG/final_scripts"

# Usage
## after setting up the variables in line 2-5, do
## chmod u+x ./Submit.sh
## ./Submit.sh




export LANG=""
cd ${workingdir}
module use /g/data/if89/apps/modulefiles
module load perllib/v5.26.3 genometools/1.6.2
module load seqtk/1.3
module load AGAT/1.4.0
module load pythonlib/3.9.2
module load cdhit/4.8.1
module load bedtools/2.31.0
mkdir -p ${workingdir}/TrainingGene
mkdir -p ${workingdir}/TrainingGene/Workingdir
mkdir -p ${workingdir}/TrainingGene/Workingdir/highQ_cDNA
mkdir -p ${workingdir}/TrainingGene/Workingdir/highQ_CDS
mkdir -p ${workingdir}/TrainingGene/Workingdir/Alignment_cDNA
mkdir -p ${workingdir}/TrainingGene/Workingdir/Alignment_CDS
mkdir -p ${workingdir}/TrainingGene/Workingdir/tmp
mkdir -p ${workingdir}/TrainingGene/Workingdir/Merged
mkdir -p ${workingdir}/Hints
mkdir -p ${workingdir}/Hints/Workingdir
mkdir -p ${workingdir}/Hints/Workingdir/weakerQ_cDNA
mkdir -p ${workingdir}/Hints/Workingdir/weakerQ_CDS
mkdir -p ${workingdir}/Hints/Workingdir/Alignment_cDNA
mkdir -p ${workingdir}/Hints/Workingdir/Alignment_CDS
mkdir -p ${workingdir}/Hints/Workingdir/tmp
mkdir -p ${workingdir}/Hints/Workingdir/Merged
mkdir -p ${workingdir}/log
mkdir -p ${workingdir}/log/PBS

grep '>' ${TRANSLATION_OUTPUT}/*.cds.all.fa | \
perl -lne '$_=~/>(\S+) len=(\d+) .* (\S+)/; $id = $1; $transcriptLen = $2; @info = split(":",$3); if ($info[3] eq "yes" and $info[4] eq "yes" and $info[5] >= 0.95 and $info[5] <= 1.05) { $id =~ /(\S+)_i\d+/; $gid=$1; print "$gid\t$id\t".join("\t",@info)."\t".($info[1]-$info[0])."\t$transcriptLen"; }' | \
sort -k1,1 -k9,9nr -k10,10nr | perl -lne '@a=split("\t",$_); unless (exists $sel{$a[0]}) { print join("\t",@a); $sel{$a[0]}="" }' > ${workingdir}/TrainingGene/Workingdir/highQ.filtered.transcripts.tsv

for i in $(ls ${TRANSLATION_OUTPUT}/*.cds.all.fa); do
    seqtk subseq -l 60 ${i} <(cut -f2 ${workingdir}/TrainingGene/Workingdir/highQ.filtered.transcripts.tsv) > ${workingdir}/TrainingGene/Workingdir/highQ_CDS/$(basename ${i} .cds.all.fa).cds.highQ.fa
done 2>&1 | tee > ${workingdir}/log/filter.log
for i in $(ls ${TRANSLATION_OUTPUT}/*.cdna.all.fa); do
    seqtk subseq -l 60 ${i} <(cut -f2 ${workingdir}/TrainingGene/Workingdir/highQ.filtered.transcripts.tsv) > ${workingdir}/TrainingGene/Workingdir/highQ_cDNA/$(basename ${i} .cdna.all.fa).cdna.highQ.fa
done 2>&1 | tee >> ${workingdir}/log/filter.log

grep '>' ${TRANSLATION_OUTPUT}/*.cds.all.fa | \
perl -lne '$_=~/>(\S+) len=(\d+) .* (\S+)/; $id = $1; $transcriptLen = $2; @info = split(":",$3); if ($info[3] eq "yes" and $info[4] eq "yes" and ($info[5] < 0.95 or $info[5] > 1.05)) { $id =~ /(\S+)_i\d+/; $gid=$1; print "$gid\t$id\t".join("\t",@info)."\t".($info[1]-$info[0])."\t$transcriptLen"; }' | \
sort -k1,1 -k9,9nr -k10,10nr | perl -lne '@a=split("\t",$_); unless (exists $sel{$a[0]}) { print join("\t",@a); $sel{$a[0]}="" }' > ${workingdir}/Hints/Workingdir/weakerQ.filtered.transcripts.tsv

for i in $(ls ${TRANSLATION_OUTPUT}/*.cds.all.fa); do
    seqtk subseq -l 60 ${i} <(cut -f2 ${workingdir}/Hints/Workingdir/weakerQ.filtered.transcripts.tsv) > ${workingdir}/Hints/Workingdir/weakerQ_CDS/$(basename ${i} .cds.all.fa).cds.weakerQ.fa
done 2>&1 | tee > ${workingdir}/log/filterWeaker.log
for i in $(ls ${TRANSLATION_OUTPUT}/*.cdna.all.fa); do
    seqtk subseq -l 60 ${i} <(cut -f2 ${workingdir}/Hints/Workingdir/weakerQ.filtered.transcripts.tsv) > ${workingdir}/Hints/Workingdir/weakerQ_cDNA/$(basename ${i} .cdna.all.fa).cdna.weakerQ.fa
done 2>&1 | tee >> ${workingdir}/log/filterWeaker.log


# Count number of jobs
export var1=$(for i in $(ls ${workingdir}/TrainingGene/Workingdir/highQ_CDS/*.cds.highQ.fa); do echo ${i}; done | wc -l)
export var2=$(for i in $(ls ${workingdir}/TrainingGene/Workingdir/highQ_cDNA/*.cdna.highQ.fa); do echo ${i}; done | wc -l)
export var3=$(for i in $(ls ${workingdir}/Hints/Workingdir/weakerQ_CDS/*.cds.weakerQ.fa); do echo ${i}; done | wc -l)
export var4=$(for i in $(ls ${workingdir}/Hints/Workingdir/weakerQ_cDNA/*.cdna.weakerQ.fa); do echo ${i}; done | wc -l)
export total=$((var1 + var2 + var3 + var4))

# Submit script that runs after minimap2 finish
export DEPEND_FOR_MINIMAP2=$(qsub -W depend=on:${total} -P ${PROJECT} -o ${workingdir}/log -v TRANSLATION_OUTPUT=${TRANSLATION_OUTPUT},workingdir=${workingdir},genome=${genome},scriptdir=${scriptdir} ${scriptdir}/1AutoSubmission/PostMinimap2.sh)

# Submit minimap2 jobs
export output=${workingdir}/TrainingGene/Workingdir/Alignment_CDS
for i in $(ls ${workingdir}/TrainingGene/Workingdir/highQ_CDS/*.cds.highQ.fa); do
    qsub -W depend=beforeok:${DEPEND_FOR_MINIMAP2} -P ${PROJECT} -o ${workingdir}/log/PBS/$(basename ${i} .fa).OU -v genome=${genome},transcriptome=${i},output=${output} ${scriptdir}/minimap2_transcriptome2genome.sh
done
export output=${workingdir}/TrainingGene/Workingdir/Alignment_cDNA
for i in $(ls ${workingdir}/TrainingGene/Workingdir/highQ_cDNA/*.cdna.highQ.fa); do
    qsub -W depend=beforeok:${DEPEND_FOR_MINIMAP2} -P ${PROJECT} -o ${workingdir}/log/PBS/$(basename ${i} .fa).OU -v genome=${genome},transcriptome=${i},output=${output} ${scriptdir}/minimap2_transcriptome2genome.sh
done
##
export output=${workingdir}/Hints/Workingdir/Alignment_CDS
for i in $(ls ${workingdir}/Hints/Workingdir/weakerQ_CDS/*.cds.weakerQ.fa); do
    qsub -W depend=beforeok:${DEPEND_FOR_MINIMAP2} -P ${PROJECT} -o ${workingdir}/log/PBS/$(basename ${i} .fa).OU -v genome=${genome},transcriptome=${i},output=${output} ${scriptdir}/minimap2_transcriptome2genome.sh
done
export output=${workingdir}/Hints/Workingdir/Alignment_cDNA
for i in $(ls ${workingdir}/Hints/Workingdir/weakerQ_cDNA/*.cdna.weakerQ.fa); do
    qsub -W depend=beforeok:${DEPEND_FOR_MINIMAP2} -P ${PROJECT} -o ${workingdir}/log/PBS/$(basename ${i} .fa).OU -v genome=${genome},transcriptome=${i},output=${output} ${scriptdir}/minimap2_transcriptome2genome.sh
done
