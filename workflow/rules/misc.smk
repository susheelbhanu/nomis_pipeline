# MAGs + gRodon + PopCOGent + antismash
# Runs miscellaneous tools for analysing growth_rate (gRodon), population diversities (PopCOGent) and biosynthetic gene clusters (BGCs; antismash)
# TODO: PopCOGent will be added in the future for "strains" of interest

localrules: install_gRodon, preprocess, merge_gRodon

###########################
# default

rule miscellaneous:
    input:
#        os.path.join(RESULTS_DIR, "bins_antismash/{sample}__{i}/knownclusterblastoutput.txt"),
        os.path.join(RESULTS_DIR, "logs/gRodon_collection.done"),
        os.path.join(RESULTS_DIR, "gRodon/merged_all_growth_prediction.txt")
    output:
        touch("status/misc.done")


####################################
# rules for Miscellaneous analyses #
####################################
# antiSMASH 
rule bins_antismash:
    input:
        os.path.join(RESULTS_DIR, "mantis_bins/{sample}__{i}.contigs.fa")
    output:
        os.path.join(RESULTS_DIR, "bins_antismash/{sample}__{i}/knownclusterblastoutput.txt")
    conda:
        os.path.join(ENV_DIR, "antismash.yaml")
    threads:
        config['antismash']['threads']
    log:
        os.path.join(RESULTS_DIR, "logs/antismash.{sample}__{i}.log")
    message:
        "Running antiSMASH for {wildcards.sample}"
    shell:
        "(date && antismash --cpus {threads} --genefinding-tool prodigal --fullhmmer --pfam2go --asf --cb-knownclusters --clusterhmmer --cf-create-clusters {input} --output-dir $(dirname {output}) && date) &> >(tee {log})"


#################
# gRodon
# Initial Setup 
rule install_gRodon:
    output:
        done=os.path.join(RESULTS_DIR, "gRodon/gRodon.installed")
    log:
        out=os.path.join(RESULTS_DIR, "logs/setup.gRodon.log")
    conda:
        os.path.join(ENV_DIR, "gRodon.yaml")
    message:
        "Setup: install R-package gRodon"
    script:
        os.path.join(SRC_DIR, "install_gRodon.R")

# Preprocessing 
rule preprocess:
    input:
        rules.prokka.output.GFF
    output:
        os.path.join(RESULTS_DIR, "prokka/{sample}__{i}_CDS_names.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/preprocess.{sample}__{i}.log")
    message:
        "Preprocessing GFFs from {wildcards.sample}"
    shell:
        """(date && sed -n '/##FASTA/q;p' {input} | awk '$3=="CDS"' | awk '{{print $9}}' | awk 'gsub(";.*","")' | awk 'gsub("ID=","")' > {output} && date) &> >(tee {log})"""

# Running gRodon 
rule gRodon:
    input:
        FFN=rules.prokka.output.FFN,
        CDS=rules.preprocess.output,
        installed=os.path.join(RESULTS_DIR, "gRodon/gRodon.installed")
    output:
        PRED=os.path.join(RESULTS_DIR, "gRodon/{sample}__{i}_growth_prediction.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/gRodon.{sample}__{i}.log")
    conda:
        os.path.join(ENV_DIR, "gRodon.yaml")
    message:
        "Growth prediction using gRodon for {wildcards.sample}"
    script:
        os.path.join(SRC_DIR, "gRodon.R")

def gRodon_aggregate(wildcards):
    checkpoint_output = checkpoints.bin_collect_mantis.get(**wildcards).output[0]
    return expand(os.path.join(RESULTS_DIR, "gRodon/{sample}__{i}_growth_prediction.txt"),
        i=glob_wildcards(os.path.join(checkpoint_output, "{i}.contigs.fa")).i
        )

rule gRodon_folder:
    input:
        gRodon_aggregate
    wildcard_constraints:
        sample="|".join(SAMPLES)
    output:
        touch(os.path.join(RESULTS_DIR, "logs/{sample}_gRodon.done"))

rule gRodon_collection_folder:
    input:
        expand(os.path.join(RESULTS_DIR, "logs/{sample}_gRodon.done"), sample=SAMPLES)
    output:
        touch(os.path.join(RESULTS_DIR, "logs/gRodon_collection.done"))


# Merging all gRodon output files
rule merge_gRodon:
    input:
        PRED=os.path.join(RESULTS_DIR, "gRodon/gRodon.installed")
    output:
        DF=os.path.join(RESULTS_DIR, "gRodon/merged_all_growth_prediction.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/gRodon.merged.log")
    conda:
        os.path.join(ENV_DIR, "r-conda.yaml")
    message:
        "Merging gRodon output for all samples"
    script:
        os.path.join(SRC_DIR, "merge_gRodon.R") 