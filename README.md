# Annotation_AusARG
Repository to keep track of changes we are making to optimise the annotation scripts for the AusARG automatic genome assembly pipeline

### Things to add/change:
- Use Diamond instead of blast to do alignment of unfiltered transcriptome to uniprot-swissprot database.
- Adjust blastxtranslation for Diamond output format.
- Adjust blastxtranslation so it prints CDS/cDNA length information in fasta header of all output file for easy UTR calculation.
- New minimap2 script to align transcriptome to genome with splice awareness. (instead of exonerate)
- New script to convert minimap2 bam format output to gff3 output.
- New script to calculate/infer UTR coordinates in genome based on cDNA alignment coordinates.
- New script to filter for high-quality transcripts based on Diamond alignment and blastxtranslation results.
- New script to get rid of redundant transcripts between different transcriptome generated from different samples. (cdhit-est)
- New script to convert gff3 output into augustus hints.
- Adjust or Make a new custom Augustus configuration file to accomodate our new hints. (extrinsic.cfg)

### Things to note
- /scripts contain scripts that are being worked on.
- /final_scripts contain scripts that are considered finished.
- /test_data contain some test data.
