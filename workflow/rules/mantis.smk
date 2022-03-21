# mantis workflow
# Runs MANTIS on the MAGs for metabolic assessment

import os
import glob
import pandas as pd

localrules: bin_collect_mantis, bin_link_mantis, mantis_config, mantis_reformat_consensus, bin_folder_sample_mantis, bin_folder_sample

###########################
# default

rule mantis:
    input:
        os.path.join(RESULTS_DIR, "logs/mantis.done")
    output:
        touch("status/mantis.done")


########################################
# checkpoint rules for collecting bins #
########################################
checkpoint bin_collect_mantis:
    input:
        os.path.join(DATA_DIR, "{sample}/run1/Binning/selected_DASTool_bins")
    output:
        directory(os.path.join(RESULTS_DIR, "mantis_links/{sample}_bins"))
    shell:
        "ln -vs {input} {output}"

rule bin_link_mantis:
    input:
        os.path.join(RESULTS_DIR, "mantis_links/{sample}_bins/{i}.contigs.fa"),
    output:
        os.path.join(RESULTS_DIR, "mantis_bins/{sample}__{i}.contigs.fa")
    wildcard_constraints:
        sample="|".join(SAMPLES)
        #i="(\w\.)+"
    shell:
        "ln -vs {input} {output}"


####################
# rules for MANTIS #
####################
# MAG list - creating a list of all the mags
# Prokka 
rule prokka:
    input:
        os.path.join(RESULTS_DIR, "mantis_bins/{sample}__{i}.contigs.fa")
    output:
        FAA=os.path.join(RESULTS_DIR, "prokka/{sample}__{i}.faa"),
        FFN=os.path.join(RESULTS_DIR, "prokka/{sample}__{i}.ffn"),
        GFF=os.path.join(RESULTS_DIR, "prokka/{sample}__{i}.gff")
    log:
        os.path.join(RESULTS_DIR, "logs/prokka.{sample}__{i}.log")
    threads:
        config['prokka']['threads']
    conda:
        os.path.join(ENV_DIR, "prokka.yaml")
    message:
        "Running Prokka on {wildcards.mag}"
    shell:
        "(date && prokka --outdir $(dirname {output.FAA}) {input} --cpus {threads} --force && date) &> >(tee {log})"

##########
# MANTIS #
rule mantis_metadata:
    input:
        txt=os.path.join(RESULTS_DIR, "data/mag_list.txt"),
        FAA=rules.prokka.output
    output:
        os.path.join(RESULTS_DIR, "data/mantis_metadata.tsv")
    shell:
        "for fname in {input.txt} ; do echo echo \"${{fname}}\"\"    \"$(echo {input.FAA}) ; done > {output}"

# Mantis: create config file from IMP config
rule mantis_config:
    output:
        "mantis.config"
    message:
        "Creating config for MANTIS"
    run:
        with open(output[0], "w") as ofile:
            # default HMMs
            for hmm_name, hmm_path in config["mantis"]["default"].items():
                ofile.write("%s=%s\n" % (hmm_name, hmm_path))
            # custom HMMs
            for hmm_path in config["mantis"]["custom"]:
                ofile.write("custom_hmm=%s\n" % hmm_path)
            # weights
            for weights_name, weights_value in config["mantis"]["weights"].items():
                ofile.write("%s=%f\n" % (weights_name, weights_value))

# Mantis: protein annotation
# NOTE: check installation before use: python submodules/mantis/ check_installation
rule mantis_run:
    input:
        FAA=rules.prokka.output,
        config="mantis.config"
    output:
        os.path.join(RESULTS_DIR, "mantis/{sample}__{i}/consensus_annotation.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/{sample}__{i}.analysis_mantis.log")
    threads:
        config['mantis']['cores']
    params:
         # The below acquires the number of cores to use from the snakemake launch command
#        CORES=int(os.environ.get("CORES"))
        cores=config['mantis']['hmmer_threads']
    conda:
        os.path.join(ENV_DIR, "IMP_MANTIS.yaml")
    message:
        "Annotation with MANTIS for {wildcards.mag}"
    shell:
        "(date && python {config[mantis][path]}/ run_mantis -t {input.FAA} --output_folder $(dirname {output}) --mantis_config {input.config} --hmmer_threads {params.cores} --cores {threads} --memory {config[mantis][single_mem]} --kegg_matrix && date) &> >(tee {log})"

# Mantis: reformat consensus output (to be imported in Python/R)
rule mantis_reformat_consensus:
    input:
        rules.mantis_run.output
    output:
        os.path.join(RESULTS_DIR, "mantis/{sample}__{i}/consensus_annotation.reformatted.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/{sample}__{i}.reformat.log")
    message:
        "Reformatting MANTIS output for {wildcards.mag}"
    shell:
        # merge annotations after "|", remove "|" separator
        "(date && paste <(cut -d '|' -f1 {input} | sed 's/\\t$//') <(cut -d '|' -f2 {input} | sed 's/^\\t//;s/\\t/;/g') > {output} && date) &> >(tee {log})"

def bins_mantis(wildcards):
    checkpoint_output = checkpoints.bin_collect_mantis.get(**wildcards).output[0]
    return expand(os.path.join(RESULTS_DIR, "mantis/{sample}__{i}/consensus_annotation.reformatted.tsv"),
        i=glob_wildcards(os.path.join(checkpoint_output, "{i}.contigs.fa")).i
        )

rule bin_folder_sample_mantis:
    input:
        bins_mantis
    wildcard_constraints:
        sample="|".join(SAMPLES)
    output:
        touch(os.path.join(RESULTS_DIR, "logs/{sample}_mantis.done"))

rule bin_folder_mantis:
    input:
        expand(os.path.join(RESULTS_DIR, "logs/{sample}_mantis.done"), sample=SAMPLES)
    output:
        touch(os.path.join(RESULTS_DIR, "logs/mantis.done"))


# # # LEGACY rules for multi-sample MANTIS, EUKARYOTES AND INTERACTIONS
#rule multi_mantis:
#    input:
#        TSV=os.path.join(RESULTS_DIR, "data/mantis_metadata.tsv"),
#        config="mantis.config"
#    output:
#        os.path.join(RESULTS_DIR, "mantis/{sample}__{i}/consensus_annotation.tsv")
#    log:
#        os.path.join(RESULTS_DIR, "logs/{sample}__{i}.analysis_mantis.log")
#    threads:
#        config['mantis']['hmmer_threads']
#    params:
#         # The below acquires the number of cores to use from the snakemake launch command
#         CORES=int(os.environ.get("CORES"))
#    conda:
#        os.path.join(ENV_DIR, "IMP_MANTIS.yaml")
#    message:
#        "Annotation with MANTIS"
#    shell:
#        "(date && python {config[mantis][path]}/ run_mantis -t {input.TSV} --output_folder $(dirname $(dirname {output})) --mantis_config {input.config} --hmmer_threads {threads} --cores {params.CORES} --memory {config[mantis][multi_mem]} --kegg_matrix && date) &> >(tee {log})"

# # rule EUK_mantis:
# #     input:
# #         FAA=os.path.join(EUK_DIR, "{euk}_good.faa"),
# #         config="mantis.config"
# #     output:
# #         os.path.join(RESULTS_DIR, "euk_mantis/{euk}/consensus_annotation.tsv")
# #     log:
# #         os.path.join(RESULTS_DIR, "logs/{euk}.analysis_mantis.log")
# #     threads:
# #         config['mantis']['hmmer_threads']
# #     params:
# #          # The below acquires the number of cores to use from the snakemake launch command
# #          CORES=int(os.environ.get("CORES"))
# #     conda:
# #         os.path.join(ENV_DIR, "IMP_MANTIS.yaml")
# #     message:
# #         "Annotation with MANTIS for {wildcards.euk}"
# #     shell:
# #         "(date && python {config[mantis][path]}/ run_mantis -t {input.FAA} --output_folder $(dirname {output}) --mantis_config {input.config} --hmmer_threads {threads} --cores {params.CORES} --memory {config[mantis][single_mem]} --kegg_matrix && date) &> >(tee {log})"

# # rule EUK_mantis_reformat_consensus:
# #     input:
# #         rules.EUK_mantis.output
# #     output:
# #         os.path.join(RESULTS_DIR, "euk_mantis/{euk}/consensus_annotation.reformatted.tsv")
# #     log:
# #         os.path.join(RESULTS_DIR, "logs/{euk}.reformat.log")
# #     message:
# #         "Reformatting MANTIS output for {wildcards.euk}"
# #     shell:
# #         # merge annotations after "|", remove "|" separator
# #         "(date && paste <(cut -d '|' -f1 {input} | sed 's/\\t$//') <(cut -d '|' -f2 {input} | sed 's/^\\t//;s/\\t/;/g') > {output} && date) &> >(tee {log})"

# # rule mantis_taxa_of_interest:
# #     input:
# #         prok=expand(os.path.join(RESULTS_DIR, "prokka/{prok}/{prok}.faa"), prok=PROK),
# #         euk=expand(os.path.join(EUK_DIR, "{euk}_good.faa"), euk=EUK),
# #         config="mantis.config"
# #     output:
# #         FAA=temp(os.path.join(RESULTS_DIR, "taxa_mantis/merged.faa")),
# #         annot=os.path.join(RESULTS_DIR, "taxa_mantis/consensus_annotation.tsv")
# #     log:
# #         os.path.join(RESULTS_DIR, "logs/taxa_of_interest_mantis.log")
# #     threads:
# #         config['mantis']['hmmer_threads']
# #     params:
# #         CORES=int(os.environ.get("CORES"))
# #     conda:
# #         os.path.join(ENV_DIR, "IMP_MANTIS.yaml")
# #     message:
# #         "Running MANTIS on taxa of interest: GL_R9_GL11_UP_2_G4.2.2 and concoct_11"
# #     shell:
# #         "(date && cat {input.prok} {input.euk} > {output.FAA} && "
# #         "python {config[mantis][path]}/ run_mantis -t {output.FAA} --output_folder $(dirname {output.annot}) --mantis_config {input.config} --hmmer_threads {threads} --cores {params.CORES} --memory {config[mantis][single_mem]} --kegg_matrix && date) &> >(tee {log})"

# # rule reformat_mantis_taxa:
# #     input:
# #         rules.mantis_taxa_of_interest.output.annot
# #     output:
# #         os.path.join(RESULTS_DIR, "taxa_mantis/consensus_annotation.reformatted.tsv")
# #     log:
# #         os.path.join(RESULTS_DIR, "logs/merged.taxa_mantis.reformat.log")
# #     message:
# #         "Reformatting the merged MANTIS output from the prok and euk bin"
# #     shell:
# #         "(date && paste <(cut -d '|' -f1 {input} | sed 's/\\t$//') <(cut -d '|' -f2 {input} | sed 's/^\\t//;s/\\t/;/g') > {output} && date) &> >(tee {log})"

