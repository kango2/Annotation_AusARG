# Pre step
## Use diamond to align unfiltered transcriptome.fa to uniprot_sprot database
## Run blastxtranslation.pl on the diamond result to get cDNA, CDS, peptide fasta file with informative headers

### This directory should contain:
### *.cds.all.fa
### *.pep.all.fa
### *.cdna.all.fa
export TRANSLATION_OUTPUT="/path/to/folder/with/blastxtranslation.pl/output"
### working directory
export workingdir="/g/data/xl04/jc4878/Annotation"
### Path to genome
export genome="/g/data/xl04/jc4878/Bassiana_publication/YAHS/genome/BASDU_HifiASM_YAHS_SUP_CurV1.1.fasta.masked"
### Path to script directory
export scriptdir="/g/data/xl04/jc4878/Annotation_AusARG/scripts"

cd ${workingdir}

### Loading modules
module use /g/data/if89/apps/modulefiles
module load perllib/v5.26.3 genometools/1.6.2
module load seqtk/1.3
module load AGAT/1.4.0
module load pythonlib/3.9.2
module load cdhit/4.8.1
module load bedtools/2.31.0

mkdir -p ${workingdir}/highQ_cDNA
mkdir -p ${workingdir}/highQ_CDS
mkdir -p ${workingdir}/log
mkdir -p ${workingdir}/log/PBS
mkdir -p ${workingdir}/log/merge
mkdir -p ${workingdir}/Alignment_cDNA
mkdir -p ${workingdir}/Alignment_CDS
mkdir -p ${workingdir}/workingdir
mkdir -p ${workingdir}/TrainingGene
mkdir -p ${workingdir}/TrainingGene/Merged
# [ Section 1 ]
## Filter for high quality transcripts, then align them to genome
### criteria:
### 1. start codon - yes
### 2. stop codon - yes
### 3. relative length to uniprot target protein : between 0.95 to 1.05
### 4. Longest isoform on cDNA level, if cDNA same length, longest on CDS level, if same CDS ength, whichever comes first
grep '>' ${TRANSLATION_OUTPUT}/*.cds.all.fa | perl -lne '$_=~/>(\S+) len=(\d+) .* (\S+)/; $id = $1; $transcriptLen = $2; @info = split(":",$3); if ($info[3] eq "yes" and $info[4] eq "yes" and $info[5] >= 0.95 and $info[5] <= 1.05) { $id =~ /(\S+)_i\d+/; $gid=$1; print "$gid\t$id\t".join("\t",@info)."\t".($info[1]-$info[0])."\t$transcriptLen"; }' | sort -k1,1 -k9,9nr -k10,10nr | perl -lne '@a=split("\t",$_); unless (exists $sel{$a[0]}) { print join("\t",@a); $sel{$a[0]}="" }' > ${workingdir}/highQ.filtered.transcripts.tsv

for i in $(ls ${TRANSLATION_OUTPUT}/*.cds.all.fa); do
    seqtk subseq -l 60 ${i} <(cut -f2 highQ.filtered.transcripts.tsv) > highQ_CDS/$(basename ${i} .cds.all.fa).cds.highQ.fa
done 2>&1 | tee > ${workingdir}/log/filter.log
for i in $(ls ${TRANSLATION_OUTPUT}/*.cdna.all.fa); do
    seqtk subseq -l 60 ${i} <(cut -f2 highQ.filtered.transcripts.tsv) > highQ_cDNA/$(basename ${i} .cdna.all.fa).cdna.highQ.fa
done 2>&1 | tee >> ${workingdir}/log/filter.log


export output=${workingdir}/Alignment_CDS
for i in $(ls ${workingdir}/highQ_CDS/*.cds.highQ.fa); do
    qsub -P ${PROJECT} -o ${workingdir}/log/PBS -v genome=${genome},transcriptome=${i},output=${output} ${scriptdir}/minimap2_transcriptome2genome.sh
done
export output=${workingdir}/Alignment_cDNA
for i in $(ls ${workingdir}/highQ_cDNA/*.cdna.highQ.fa); do
    qsub -P ${PROJECT} -o ${workingdir}/log/PBS -v genome=${genome},transcriptome=${i},output=${output} ${scriptdir}/minimap2_transcriptome2genome.sh
done

## Run below after the minimap2 job finishes, each transcriptome should only take ~5 min
for i in $(ls ${workingdir}/Alignment_cDNA/*.bam); do
  agat_convert_minimap2_bam2gff.pl -i ${i} -o ${workingdir}/Alignment_cDNA/$(basename ${i} .bam).gff3
done 2>&1 | tee > ${workingdir}/log/convertBAMtoGFF3.log
for i in $(ls ${workingdir}/Alignment_CDS/*.bam); do
  agat_convert_minimap2_bam2gff.pl -i ${i} -o ${workingdir}/Alignment_CDS/$(basename ${i} .bam).gff3
done 2>&1 | tee >> ${workingdir}/log/convertBAMtoGFF3.log

for i in $(ls ${workingdir}/Alignment_cDNA/*_2genome.gff3); do
  python ${scriptdir}/processminimap2gff3.py ${i} | \
  # Pipe output to awk command to replace ID= or Parent= with the string in column 10, get rid of column 10, then replace cDNA_match_part with exon
  awk 'BEGIN {FS="\t"; OFS="\t"} {
    if ($0 ~ /^#/) {
        print
        next
    }
    if ($3 == "cDNA_match") {
        $3="gene"
        split($9, parts, ";")
        for (i in parts) {
            split(parts[i], kv, "=")
            if (kv[1] == "ID") {
                parts[i] = kv[1]"=" $10".gene"
            }
        }
        $9 = parts[1]
        sub(/;$/, "", $9)
        print $1"\t"$2"\tmRNA\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tID="$10";Parent="$10".gene"
    } else if ($3 == "cDNA_match_part") {
        $3="exon"
        split($9, parts, ";")
        for (i in parts) {
            split(parts[i], kv, "=")
            if (kv[1] == "Parent") {
                parts[i] = kv[1]"=" $10
            }
            else if (kv[1] == "ID") {
                parts[i] = kv[1]"=" $10".exonOnLine"NR
            }
        }
        $9 = parts[1] ";" parts[2]
        sub(/;$/, "", $9)
    }
    print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9
}' > ${workingdir}/Alignment_cDNA/$(basename ${i} .gff3)_processed.gff3
done 2>&1 | tee > ${workingdir}/log/processcDNAGFF3.log
for i in $(ls ${workingdir}/Alignment_CDS/*_2genome.gff3); do
  python ${scriptdir}/processminimap2gff3.py ${i} | \
  # Pipe output to awk command to replace ID= or Parent= with the string in column 10, get rid of column 10, then replace cDNA_match_part with exon
  awk 'BEGIN {FS="\t"; OFS="\t"} {
    if ($0 ~ /^#/) {
        print
        next
    }
    if ($3 == "cDNA_match") {
        $3="gene"
        split($9, parts, ";")
        for (i in parts) {
            split(parts[i], kv, "=")
            if (kv[1] == "ID") {
                parts[i] = kv[1]"=" $10".gene"
            }
        }
        $9 = parts[1]
        sub(/;$/, "", $9)
        print $1"\t"$2"\tmRNA\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tID="$10";Parent="$10".gene"
    } else if ($3 == "cDNA_match_part") {
        $3="CDS"
        split($9, parts, ";")
        for (i in parts) {
            split(parts[i], kv, "=")
            if (kv[1] == "Parent") {
                parts[i] = kv[1]"=" $10
            }
            else if (kv[1] == "ID") {
                parts[i] = kv[1]"=" $10".CDSOnLine"NR
            }
        }
        $9 = parts[1] ";" parts[2]
        sub(/;$/, "", $9)
    }
    print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9
}' > ${workingdir}/Alignment_CDS/$(basename ${i} .gff3)_processed.gff3
done 2>&1 | tee > ${workingdir}/log/processCDSGFF3.log

# [ Section 2]
## Find transcript where the cDNA alignment completely cover 100% of the CDS alignment in bedtools intersect
## Calculate UTR from exon/CDS difference
## Infer introns for hints
for i in $(ls ${workingdir}/Alignment_cDNA/*_processed.gff3); do
    awk '$3 == "gene"' ${i} | awk 'BEGIN {FS="\t"; OFS="\t"} 
    {
        if ($3 == "gene") {
            split($9, parts, ";")
            for (i in parts) {
                split(parts[i], kv, "=")
            }
            $9 = parts[1]
            sub(/\.gene$/, "", $9)
            sub(/^ID=/, "", $9)
        }
        print $1"\t"$4"\t"$5"\t"$9"\t"$6"\t"$7
    }' > ${workingdir}/workingdir/cDNA.bed
    awk '$3 == "gene"' ${workingdir}/Alignment_CDS/$(basename ${i} cdna.highQ_2genome_processed.gff3)cds.highQ_2genome_processed.gff3 | awk 'BEGIN {FS="\t"; OFS="\t"} 
    {
        if ($3 == "gene") {
            split($9, parts, ";")
            for (i in parts) {
                split(parts[i], kv, "=")
            }
            $9 = parts[1]
            sub(/\.gene$/, "", $9)
            sub(/^ID=/, "", $9)
        }
        print $1"\t"$4"\t"$5"\t"$9"\t"$6"\t"$7
    }' > ${workingdir}/workingdir/CDS.bed
    bedtools intersect -F 1.0 -a ${workingdir}/workingdir/cDNA.bed -b ${workingdir}/workingdir/CDS.bed -wa -wb -s | awk '$4 == $10 {print $4}' > ${workingdir}/workingdir/both.txt
    echo '##gff-version 3' > ${workingdir}/workingdir/merge.gff3
    grep -w -f ${workingdir}/workingdir/both.txt ${i} >> ${workingdir}/workingdir/merge.gff3
    awk '$3 == "CDS"' ${workingdir}/Alignment_CDS/$(basename ${i} cdna.highQ_2genome_processed.gff3)cds.highQ_2genome_processed.gff3 | grep -w -f ${workingdir}/workingdir/both.txt >> ${workingdir}/workingdir/merge.gff3
    gt gff3 -sortlines -tidy -retainids ${workingdir}/workingdir/merge.gff3 | gt gff3 -retainids -sort > ${workingdir}/workingdir/merge2.gff3
    python3 ${scriptdir}/addUTRs.py ${workingdir}/workingdir/merge2.gff3 ${workingdir}/workingdir/merge2_with_UTRs.gff3
    gt gff3 -retainids -addintrons -sort ${workingdir}/workingdir/merge2_with_UTRs.gff3 | \
    awk 'BEGIN {FS="\t"; OFS="\t"} {
        if ($3 == "intron") {
            $2="genometools"
            split($9, parts, ";")
            for (i in parts) {
                split(parts[i], kv, "=")
            $9 = parts[1]".intronOn"NR";"parts[1]
            sub(/^Parent=/, "ID=", $9)
            print
            }
        }
        else {
            print
        }
    }' > ${workingdir}/workingdir/merge3_with_UTRs.gff3
    mv ${workingdir}/workingdir/merge3_with_UTRs.gff3 ${workingdir}/TrainingGene/Merged/$(basename ${i} cdna.highQ_2genome_processed.gff3)highQ_with_UTRs.gff3
done 2>&1 | tee > ${workingdir}/log/mergeGFF3.log
awk '$1 != "warning:"' ${workingdir}/log/mergeGFF3.log > ${workingdir}/log/mergeGFF3.log2

rm -rf ${workingdir}/workingdir


# [ Section 3 ]
## Getting rid of redundant protein sequences in different transcriptome using cd-hit, protein with sequence identity >= 80% are clustered and a representative sequence is given
## Note: non-redundancy is important in optimising augustus, the criteria given in the link below is >= 80% identical on aa level
### https://vcru.wisc.edu/simonlab/bioinformatics/programs/augustus/docs/tutorial2015/training.html
mkdir -p ${workingdir}/workingdir
for i in $(ls ${TRANSLATION_OUTPUT}/*.pep.all.fa); do
    seqtk subseq -l 60 ${i} <(cut -f2 ${workingdir}/highQ.filtered.transcripts.tsv) >> ${workingdir}/workingdir/transcriptome.pep.highQ.fa
done

### Change -T to ${PBS_NCPUS} if running on PBS
cd-hit -i ${workingdir}/workingdir/transcriptome.pep.highQ.fa -o ${workingdir}/workingdir/transcriptome.pep.highQ.representative.fa -c 0.8 -aS 0.9 -g 1 -d 0 -n 3 -M 75000 -T 1

seqtk seq -C ${workingdir}/workingdir/transcriptome.pep.highQ.representative.fa | grep ">" | sed 's/>//g' | sort > ${workingdir}/workingdir/transcriptome.pep.highQ.representative.list
### grep these out from gff3
grep --no-filename -w -f ${workingdir}/workingdir/transcriptome.pep.highQ.representative.list ${workingdir}/TrainingGene/Merged/*highQ_with_UTRs.gff3 > ${workingdir}/TrainingGene/Training_gene.gff3

### Find overlapping genes, get rid of one of them. 1kbp are added to either side of the gene so the 2kp flanking region generated later on will have no coding region
awk '$3 == "gene"' ${workingdir}/TrainingGene/Training_gene.gff3 | awk 'BEGIN {FS="\t"; OFS="\t"} 
{
    if ($3 == "gene") {
        split($9, parts, ";")
        for (i in parts) {
            split(parts[i], kv, "=")
        }
        $9 = parts[1]
        sub(/\.gene$/, "", $9)
        sub(/^ID=/, "", $9)
    }
    print $1"\t"$4-1000"\t"$5+1000"\t"$9"\t"$6"\t"$7
}' | awk 'BEGIN {FS="\t"; OFS="\t"} {if ($2 <= 0) {$2 = "1"} {print}}' > ${workingdir}/workingdir/Training_gene.bed

bedtools intersect -a ${workingdir}/workingdir/Training_gene.bed -b ${workingdir}/workingdir/Training_gene.bed -wa -wb | \
awk '$4 < $10 {print $4} $4 > $10 {print $10}' | sort | uniq > ${workingdir}/workingdir/overlapping.list
grep -w -v -f ${workingdir}/workingdir/overlapping.list ${workingdir}/TrainingGene/Training_gene.gff3 > ${workingdir}/TrainingGene/Training_gene_NoOverlap.gff3

### Now the same thing for hints, but without the 1kbp, also process UTR into UTRpart and exon into exonpart
awk '$3 == "gene"' ${workingdir}/TrainingGene/Training_gene.gff3 | awk 'BEGIN {FS="\t"; OFS="\t"} 
{
    if ($3 == "gene") {
        split($9, parts, ";")
        for (i in parts) {
            split(parts[i], kv, "=")
        }
        $9 = parts[1]
        sub(/\.gene$/, "", $9)
        sub(/^ID=/, "", $9)
    }
    print $1"\t"$4"\t"$5"\t"$9"\t"$6"\t"$7
}' > ${workingdir}/workingdir/Training_gene_No1Kbp.bed
bedtools intersect -a ${workingdir}/workingdir/Training_gene_No1Kbp.bed -b ${workingdir}/workingdir/Training_gene_No1Kbp.bed -wa -wb | \
awk '$4 < $10 {print $4} $4 > $10 {print $10}' | sort | uniq > ${workingdir}/workingdir/overlapping_hints.list
grep -w -v -f ${workingdir}/workingdir/overlapping_hints.list ${workingdir}/TrainingGene/Training_gene.gff3 | \
awk 'BEGIN {FS="\t"; OFS="\t"} {
        if ($3 == "intron") {
            split($9, parts, ";")
            for (i in parts) {
                split(parts[i], kv, "=")
            }
            $9 = parts[2]";pri=4;src=E"
            sub(/^Parent=/, "grp=", $9)
            print
        }
        else if ($3 == "exon") {
            $3="exonpart"
            split($9, parts, ";")
            for (i in parts) {
                split(parts[i], kv, "=")
            }
            $9 = parts[2]";pri=4;src=E"
            sub(/^Parent=/, "grp=", $9)
            print
            }
        else if ($3 == "five_prime_UTR" || $3 == "three_prime_UTR") {
            $3="UTRpart"
            split($9, parts, ";")
            for (i in parts) {
                split(parts[i], kv, "=")
            }
            $9 = parts[2]";pri=4;src=E"
            sub(/^Parent=/, "grp=", $9)
            print
        }
    }' > ${workingdir}/TrainingGene/Hints.gff3


# Continue with TrainingAugustus.sh
## Shared variables are:
### ${workingdir}
### ${genome}