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
- Adjust bonus/malus in Augustus hints configuration file to accomodate our new hints. (extrinsic.cfg)

### Things to note
- /scripts contain scripts that are being worked on.
- /final_scripts contain scripts that are considered finished.
- /test_data contain some test data.
- /extrinsic contain the default extrinsic directory content in Augustus
- /workingdir directory where I test scripts
- extrinsic.MPE.cfg is the hints bonus/malus config file I have been using, I made a copy from /extrinsic to the master directory to work on it.
- Need to adjust this config file so it make more sense, for example by default the feature "start" gets a malus of 0.8 if Augustus predicts a "start" at a region with no "start" hints supporting it. Our new hints file (or any old ones) by nature has no "start" feature in it to begin with, so all start predicted by Augustus automatically gets a 0.8 malus.
- Our hints file only has exonpart and intronpart, so these are the only two feature we should configure
- Or maybe we could change hints file to exon and intron