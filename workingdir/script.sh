# generate some test data - transcriptome (cDNA)
cd /workingdir


# /g/data/xl04/hrp561/basdurnaseq/AA045725_Brain_ZWf_renamed.fasta
# AA045725_Brain_ZWf_DN20936_c0_g2_i1
# /g/data/xl04/hrp561/basdurnaseq/mixedID_Ovary_ZWf_renamed.fasta
# mixedID_Ovary_ZWf_DN11255_c0_g2_i1

module load seqtk
seqtk subseq -l 60 /g/data/xl04/hrp561/basdurnaseq/AA045725_Brain_ZWf_renamed.fasta <(echo -e "AA045725_Brain_ZWf_DN20936_c0_g2_i1") > input/transcriptome_subset.fasta
seqtk subseq -l 60 /g/data/xl04/hrp561/basdurnaseq/mixedID_Ovary_ZWf_renamed.fasta <(echo -e "mixedID_Ovary_ZWf_DN11255_c0_g2_i1") >> input/transcriptome_subset.fasta

module unload minimap2 samtools
module load minimap2/2.26 samtools/1.19.2

export transcriptome="/g/data/xl04/jc4878/Annotation_AusARG/workingdir/input/transcriptome_subset.fasta"
export input="/g/data/xl04/jc4878/Annotation_AusARG/workingdir/input"
export output="/g/data/xl04/jc4878/Annotation_AusARG/workingdir/output"
export genome="/g/data/xl04/jc4878/Bassiana_publication/YAHS/genome/BASDU_HifiASM_YAHS_SUP_CurV1.1.fasta.masked"
export PBS_NCPUS="1"
# Check if transcriptome ends with .fa or .fasta
if [[ $transcriptome == *.fa ]]; then
    export base=$(basename "$transcriptome" .fa)
elif [[ $transcriptome == *.fasta ]]; then
    export base=$(basename "$transcriptome" .fasta)
else
    echo "Unsupported transcriptome file extension."
    exit 1
fi

minimap2 -t ${PBS_NCPUS} -ax splice:hq \
${genome} \
${transcriptome} | samtools sort -@ ${PBS_NCPUS} -O BAM -o ${output}/${base}_2genome.bam



module unload AGAT samtools perllib parallel
module load AGAT/1.4.0 parallel/20191022

export input="/g/data/xl04/jc4878/Annotation_AusARG/workingdir/output"
export output="/g/data/xl04/jc4878/Annotation_AusARG/workingdir/output"
# bam to gff3 conversion using AGAT
ls ${input}/*.bam | parallel --jobs ${PBS_NCPUS} '
  export base=$(basename {} .bam)
  agat_convert_minimap2_bam2gff.pl -i {} -o ${output}/${base}.gff3
'



# gff3 to hints.gff3
module load pythonlib/3.9.2
module load genometools/1.6.2
## We can either use exonpart/intronpart or exon/intron if we are confident about the transcripts and their alignments (and give augustus a big bonus for predicting something that follows the hint)
export input="/g/data/xl04/jc4878/Annotation_AusARG/workingdir/output/transcriptome_subset_2genome.gff3"
export base=$(basename ${input} .gff3)
export output="/g/data/xl04/jc4878/Annotation_AusARG/workingdir/output"

python /g/data/xl04/jc4878/Annotation_AusARG/scripts/processminimap2gff3.py ${input} | \
# Pipe output to awk command to replace ID= or Parent= with the string in column 10, get rid of column 10, then replace cDNA_match_part with exon
awk 'BEGIN {FS="\t"; OFS="\t"} {
    if ($3 == "cDNA_match") {
        split($9, parts, ";")
        for (i in parts) {
            split(parts[i], kv, "=")
            if (kv[1] == "ID") {
                parts[i] = kv[1]"=" $10
            }
        }
        $9 = ""
        for (i in parts) {
            $9 = $9 parts[i] ";"
        }
        sub(/;$/, "", $9)
    } else if ($3 == "cDNA_match_part") {
        split($9, parts, ";")
        for (i in parts) {
            split(parts[i], kv, "=")
            if (kv[1] == "Parent") {
                parts[i] = kv[1]"=" $10
            }
        }
        $9 = ""
        for (i in parts) {
            $9 = $9 parts[i] ";"
        }
        sub(/;$/, "", $9)
    }
    print $0
}' | cut --complement -f10 | awk 'BEGIN{FS=OFS="\t"} $3=="cDNA_match_part" {$3="exon"} 1' | \
# Pipe output to genometools to infer introns between exon features, get rid of the region lines
gt gff3 --addintrons --retainids --checkids -setsource minimap2 -addids no | grep -v "###" | \
# Pipe output to awk command, extracting only exon and intron features to modify 9th column to hint format, including grp=, pri4, and src=E
awk -F'\t' -v OFS='\t' '$3 == "exon" || $3 == "intron" { match($9, /Parent=([^;\t ]+)/, arr); $9 = "grp=" arr[1] ";pri=4;src=E"; print }' | \
# Pipe output to awk command to replace exon with exonpart, intron with intronpart
awk 'BEGIN{FS=OFS="\t"} $3=="exon" {$3="exonpart"} 1' | awk 'BEGIN{FS=OFS="\t"} $3=="intron" {$3="intronpart"} 1' > ${output}/${base}.hints.gff3

