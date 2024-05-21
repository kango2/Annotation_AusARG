# Annotation_AusARG

Gene annotation pipeline for AusARG genomes. Scripts here are designed to annotate new genomes generated from the AusARG assembly pipeline. The annotation pipeline is divided into three main parts/scripts:
1. **GenerateTrainingGene** Generate high quality training genes for training Augustus, and generate Augustus hints.
2. **TrainingAugustus** Generate new Augustus species-specific parameters and optimise it
3. **RunningAugustus** Run gene prediction using trained parameters

# Usage

To run the pipeline for the first time on a new species, follow this progression:  
1 -> 2 -> 3  
To run the pipeline for a different genome of the same species (after initial run), follow this progression:  
1 -> 3  

# 1 Generate Training Gene
~30 minutes to run  
Best to run it in interactive job or manually set $LANG to "" or the result could be different due to sorting  

**Required Inputs (more info in script file):**
1. path to blastxtranslation outputs
2. path to working directory
3. path to soft-masked genome
4. path to script directory

**Output:**
```
├── ${workingdir}
│   └── TrainingGene
│       └── Training_gene_NoOverlap.gff3
│   └── Hints
│       └── Hints.gff3
```

# 2 Training Augustus
~30 hours to run  

**Required Inputs (more info in script file):**
1. name of species (with no whitespace)
2. path to working directory
3. path to soft-masked genome

**Output:**
```
├── ${workingdir}
│   └── Augustus
│       └── config (contains trained parameter for the new species)
```

# 3 Running Augustus
~4 hours to run  

**Required Inputs (more info in script file):**
1. name of species (with no whitespace)
2. path to working directory
3. path to soft-masked genome
4. path to Augustus config (Except initial run)
5. path to custom extrinsic file
6. path to uniprot-swissprot diamond database
7. path to uniprot-trembl diamond database

**Output:**
```
├── ${workingdir}/Augustus_annotation.gff3
├── ${workingdir}/Augustus_peptide.fasta
├── ${workingdir}/Augustus_gene_table.tabular (in TSV format)
```