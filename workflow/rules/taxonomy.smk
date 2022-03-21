# taxonomic workflow
# Runs GTDBtk & CheckM on the bins

###########################
# default

rule taxonomy:
    input:
        os.path.join(RESULTS_DIR, "gtdbtk_output"),
        os.path.join(RESULTS_DIR, "checkm_output")
    output:
        touch("status/gtdbtk_checkm.done")


########################################
# rules for Taxonomy and completenesss #
########################################
# GTDBTK taxonomy
rule gtdbtk:
    input:
        os.path.join(RESULTS_DIR, "bins/bin_collection.done")
    output:
        directory(os.path.join(RESULTS_DIR, "gtdbtk_output"))
    log:
        os.path.join(RESULTS_DIR, "logs/gtdbtk.log")
    conda:
        os.path.join(ENV_DIR, "gtdbtk.yaml")
    params:
        config["gtdbtk"]["path"]
    threads:
        config["gtdbtk"]["threads"]
    message:
        "Running GTDB toolkit on MAGs"
    shell:
        "(date && export GTDBTK_DATA_PATH={params} && gtdbtk classify_wf --cpus {threads} -x fa --genome_dir $(dirname {input}) --out_dir {output} && date) &> >(tee {log})"

# Checking bin quality
rule checkm:
    input:
        os.path.join(RESULTS_DIR, "bins/bin_collection.done")
    output:
        directory(os.path.join(RESULTS_DIR, "checkm_output"))
    log:
        os.path.join(RESULTS_DIR, "logs/checkm.log")
    conda:
        os.path.join(ENV_DIR, "checkm.yaml")
    threads:
        config["checkm"]["threads"]
    message:
        "Checking MAG quality"
    shell:
        "(date && checkm lineage_wf -r -t {threads} -x fa $(dirname {input}) {output} && date) &> >(tee {log})"