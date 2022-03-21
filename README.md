# nomis_pipeline
## About
- Repository containing workflows for `IMP3` downstream analyses
- Related project(s): [NOMIS](https://nomis-data.epfl.ch/)

# Setup

## Conda

[Conda user guide](https://docs.conda.io/projects/conda/en/latest/user-guide/index.html)

```bash
# install miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod u+x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh # follow the instructions
```

Getting the repository including sub-modules
```bash
git clone --recurse-submodules ssh://git@git-r3lab-server.uni.lu:8022/susheel.busi/nomis_pipeline.git
```

Create the main `snakemake` environment

```bash
# create venv
conda env create -f requirements.yml -n "snakemake"
```

## Dependencies

The successful completion requires tools created by others

1. [MANTIS](https://github.com/PedroMTQ/mantis)
2. [MAGICCAVE/MAGICLAMP](https://github.com/Arkadiy-Garber/MagicLamp)
  
Notes: 
- Dependencies are included as `submodules` where possible
- However, installation issues may persist
- If so, check the respective repositories listed
  
## How to run

The workflow can be launched using one of the option as follows

```bash
./config/sbatch.sh
```

(or)

```bash
CORES=48 snakemake -s workflow/Snakefile --configfile config/config.yaml --use-conda --conda-prefix ${CONDA_PREFIX}/pipeline --cores $CORES -rpn
```
(or)

Note: For running on `esb-compute-01` or `litcrit`  adjust the `CORES` as needed to prevent `MANTIS` from spawning too many workers and launch as below

```bash
CORES=24 snakemake -s workflow/Snakefile --configfile config/config.yaml --use-conda --conda-prefix ${CONDA_PREFIX}/pipeline --cores $CORES -rpn
```

## Configs

All config files are stored in the folder `config/`:

# Workflows

1. `imp workflow`: setup the folders required for running IMP3 on each sample
2. `viruses workflow`: run VIBRANT and vCONTACT2 on assemblies, including CheckV on vibrant_output
3. `eukaryotes workflow`: running EUKUlele on assemblies 
4. `bins workflow`: collects all bins together for taxonomy analyses
5. `taxonomy workflow`: run GTDBtk and CheckM on the bins
6. `functions workflow`: runs METABOLIC, MAGICCAVE and FUNCS analyses
7. `mantis workflow`: runs MANTIS on the bins explicitly
8. `euk_bin workflow`: performs coassembly specifically for eukaryotes (EukRep) and runs binning with CONCOCT
9. `coassembly_binning`: performs coassembly for all samples and subsequent binning
10. `misc workflow`: runs gRodon, antismash. To be implemented PopCOGent,and potentially anvi'o coassembly/binning.

Relevant paremters which have to be changed are listed for each workflow and config file.
Parameters defining system-relevant settings are not listed but should be also be changed if required, e.g. number of threads used by certain tools etc.

## STEPS

The workflow is setup in multiple steps. Prior to running change the following

- config: `config/`
  - `config.yaml`:
    - change `steps`

Options:

1. *imp*
2. *viruses*
3. *eukaryotes*
4. *bins*
5. *taxonomy*
6. *functions*
7. *mantis*
8. *euk_bin*
9. *coassembly_binning*
10. *misc*

IMPORTANT NOTE: only the `imp` step should be run first, followed by `launching IMP3` outside of this pipeline. Subsequent other `STEPS` can be run

### Launching `IMP3`

Per-sample IMP3 can be launched as follows:

```bash
chmod -R 775 ${SAMPLE}      # adding permissions
cd ${SAMPLE}
sbatch ./launchIMP.sh       # on IRIS
```


## imp workflow

Download raw data required for the analysis.

- config: `config/`
  - `config.yaml`:
    - change `work_dir`
  - `sbatch.sh`
    - change `SMK_ENV`
    - if not using `slurm` to submit jobs remove `--cluster-config`, `--cluster` from the `snakemake` CMD
  - `slurm.yaml` (only relevant if using `slurm` for job submission)
- workflow: `workflow/`

Prior to running the `imp workflow` make the following adjustments.

- IMP_config.yaml: 
  - `workflow/notes/IMP_config.yaml`
    - change `Metagenomics`
- run_threads:
  - `workflow/notes/runIMP.sh`
    - change: `threads`
- launch_threads:
  - `workflow/notes/runIMP.sh`
    - change: `-n8`

IMPORTANT Note:

This above workflow should be run first, followed by launching IMP3 outside of this pipeline and then subsequent `STEPS` can be run

## Main workflow

Main analysis workflow: given SR FASTQ files, run all the steps to generate required output. This includes:

- setting up folders for IMP
- viral and eukaryotic annotations
- functional analyses and
- taxonomic analyses (optional)

The workflow is run per sample and might require a couple of days to run depending on the sample, used configuration and available computational resources.
Note that the workflow will create additional output files not necessarily required to re-create the figures shown in the manuscript.

- config:
  - per sample
  - `config/<sample>/config.yaml`
    - change all path parameters (not all databases are required, see above)
  - `config/<sample>/sbatch.yaml`
    - change `SMK_ENV`
    - if not using `slurm` to submit jobs remove `--cluster-config`, `--cluster` from the `snakemake` CMD
  - `config/<sample>/slurm.yaml` (only relevant if using `slurm` for job submission)
- workflow: `workflow/`

## Report workflow (2021-05-26 15:54:59: NOT implemented)

This workflow creates various summary files, plots and an HTML report for a sample using the output of the main workflow.

Note: How the metaP peptide/protein reports were generated from raw metaP data is described in `notes/gdb_metap.md`.

- config:
  - sample configs used for the main workflow
- workflow: `workflow_report/`

To execute this workflow for all samples:
```bash
./config/reports.sh "YourEnvName" "WhereToCreateCondEnvs"
```

## Figures workflow (2021-05-26 15:54:52: NOT implemented)

Re-create figures (and tables) used in the manuscript.
This workflow should be only run after running the main workflow and report workflow for all samples.

- config: `config/fig.yaml`
  - change `work_dir`
  - change paths for all samples in `samples`
- workflow: `workflow_figures/`

```bash
conda activate "YourEnvName"
snakemake -s workflow_figures/Snakefile --cores 1 --configfile config/fig.yaml --use-conda --conda-prefix "WhereToCreateCondEnvs" -rpn # dry-run
```

# Notes

Notes for manual/additional analyses done using the generated data.
