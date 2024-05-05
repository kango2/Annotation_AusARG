#!/bin/bash
#PBS -N Augustus_training
#PBS -l ncpus=16,walltime=48:00:00,storage=gdata/if89+gdata/xl04,mem=50GB,jobfs=50GB
#PBS -j oe


### Loading modules
module use /g/data/if89/apps/modulefiles
module load Augustus/3.4.0 perllib/v5.26.3 blat/37 RepeatMasker/4.1.2-p1 scipio/1.4 pblat/2.5 pslCDnaFilter/0 parallel/20191022 exonerate/2.2.0

### Copying augustus config files to a writable space and setting up path
cd ${workingdir}
mkdir -p ${workingdir}/Augustus
mkdir -p ${workingdir}/Augustus/training
mkdir -p ${workingdir}/Augustus/config


rsync -a $AUGUSTUS_CONFIG_PATH/ Augustus/config/
export AUGUSTUS_CONFIG_PATH=${workingdir}/Augustus/config

# ${training} is a new variable that should be the path to the 500 manually selected training genes in gff3 format, containing the features CDS, three_prime_utr and five_prime_utr in 3rd column
cd ${workingdir}/Augustus/training
gff2gbSmallDNA.pl --softmasked ${training} ${genome} 2000 ${workingdir}/Augustus/training/training.gb

# Split the 500 into 200 (.test) for evaluation test and 300 (.train) for training
randomSplit.pl training.gb 200
# into training.gb.test (200)
# and training.gb.train (300)

# First we make parameters files for the new species
new_species.pl --species=${species}
# Concatenate generic metaparameters with UTR metaparameters
grep -h -vE '^(#|$)' ${workingdir}/Augustus/config/species/${species}/${species}_metapars.cfg \
${workingdir}/Augustus/config/species/${species}/${species}_metapars.utr.cfg \
> ${workingdir}/Augustus/config/species/${species}/${species}_metapars_and_utr.cfg

# an initial training using genes not selected for optimising, this is fast
etraining --species=${species} --UTR=on --print_utr=on --stopCodonExcludedFromCDS=f ${workingdir}/Augustus/training/training.gb.train
augustus --species=${species} --UTR=on --print_utr=on --stopCodonExcludedFromCDS=f ${workingdir}/Augustus/training/training.gb.test | tee first_evaluation.out
grep -A 36 Evaluation first_evaluation.out > ${workingdir}/first_evaluation.report

# Now we optimize with 500 genes, 200 for evaluation and all 500 for training, the max 5 rounds has been chosen but it will finish earlier if no improvement are found
optimize_augustus.pl \
--species=${species} \
--cpus=${PBS_NCPUS} \
--rounds=5 \
${workingdir}/Augustus/training/training.gb.test \
--onlytrain=${workingdir}/Augustus/training/training.gb.train \
--UTR=on \
--metapars=${workingdir}/Augustus/config/species/${species}/${species}_metapars_and_utr.cfg \
--cleanup=1 \
> ${workingdir}/OptimiseAugustus.running

# Retrain after optimization
etraining --species=${species} --UTR=on --print_utr=on --stopCodonExcludedFromCDS=f ${workingdir}/Augustus/training/training.gb.train

# final evaluation
augustus --species=${species} --UTR=on --print_utr=on --stopCodonExcludedFromCDS=f ${workingdir}/Augustus/training/training.gb.test | tee final_evaluation.out
grep -A 36 final_evaluation.out > ${workingdir}/final_evaluation.report



# rename log file from .running to .done
mv ${workingdir}/OptimiseAugustus.running ${workingdir}/OptimiseAugustus.done
