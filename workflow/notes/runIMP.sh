#!/bin/bash -l

export PATH=$PATH:/work/projects/ecosystem_biology/local_tools/bin
export PATH=/work/projects/ecosystem_biology/local_tools/minicondaESB/condabin:/work/projects/ecosystem_biology/local_tools/bin:$PATH

IMPROOT=/work/projects/ecosystem_biology/local_tools/IMP3
THREADS=8
TMPDIR=tmp

snakemake -s $IMPROOT/Snakefile --configfile "$sample"_config.yaml -j $THREADS --use-conda --conda-prefix $IMPROOT/conda --unlock --rerun-incomplete
snakemake -s $IMPROOT/Snakefile --configfile "$sample"_config.yaml -j $THREADS --use-conda --conda-prefix $IMPROOT/conda --rerun-incomplete

# moving the files to a different location
# rsync --remove-source-files -auv --no-p --no-g -P /scratch/users/sbusi/metaG_JULY_2020/IMP3/"$sample" /work/projects/nomis/metaG_JULY_2020/IMP3/.

### NOTE to self ###
# Use the below to check for 'all.done' file and then move if exists
# [ -f "$sample"/run1/status/all.done ] && { echo ""$sample" is IMPed"; rsync --remove-source-files -auv --no-p --no-g -P /scratch/users/sbusi/metaG_JULY_2020/IMP3/"$sample" /work/projects/nomis/metaG_JULY_2020/IMP3/.; }
