# euk_bin workflow
# Runs several steps to obtain Eukaryotic bins

###########################
# default
rule euk_bin: 
    input: 
        os.path.join(RESULTS_DIR, "euk_bin/concoct/bins")
    output:
        touch("status/euk_bin.done")


#################
# Preprocessing #
#################
rule deduplicate:
    input:
        r1=os.path.join(DATA_DIR, "{sample}/run1/Preprocessing/mg.r1.preprocessed.fq"),
        r2=os.path.join(DATA_DIR, "{sample}/run1/Preprocessing/mg.r2.preprocessed.fq")
    output:
        odup1=os.path.join(RESULTS_DIR, "euk_bin/dedup/{sample}_R1.fastq.gz"), 
        odup2=os.path.join(RESULTS_DIR, "euk_bin/dedup/{sample}_R2.fastq.gz")
    log:
        os.path.join(RESULTS_DIR, "logs/dedup.{sample}.log")
    threads:
        config["clumpify"]["threads"]
    conda:
        os.path.join(ENV_DIR, "bbmap.yaml")
    message:
        "Removing duplicate reads for easier downstream assembly"
    shell:
        "(date && clumpify.sh in={input.r1} in2={input.r2} out={output.odup1} out2={output.odup2} dupedist={config[clumpify][dupedist]} dedupe=t optical=t threads={threads} groups={config[clumpify][groups]} -Xmx{config[clumpify][memory]} && date) &> >(tee {log})"

###################
# Read extraction #
###################
# Extracting only those reads classified as Eukaryotes
# Taxonomy with Kraken2
rule euk_bin_kraken2:
    input:
        dedup1=rules.deduplicate.output.odup1,
        dedup2=rules.deduplicate.output.odup2,
        database=config['kraken2']['db']
    output:
        rep=os.path.join(RESULTS_DIR, "euk_bin/kraken2/{sample}.kraken.report.txt"),
        summary=os.path.join(RESULTS_DIR, "euk_bin/kraken2/{sample}.kraken.summary.out")
    log:
        os.path.join(RESULTS_DIR, "logs/kraken2.euk_bin.{sample}.log")
    threads:
        config['kraken2']['threads']
    conda:
        os.path.join(ENV_DIR, "kraken2.yaml")
    message:
        "Kraken2 taxonomy for eukaryotic binning"
    shell:
        "(date && kraken2 --threads {threads} --db {input.database} --use-names --confidence 0.5 --paired {input.dedup1} {input.dedup2} --gzip-compressed --output {output.summary} --report {output.rep} && date) &> >(tee {log})"

rule headers:
    input:
        rules.euk_bin_kraken2.output.summary
    output:
        os.path.join(RESULTS_DIR, "euk_bin/headers/{sample}.headers.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/headers.euk_bin.{sample}.log")
    message:
        "Extracting 'unclassified' read headers"
    shell:
        "(date && awk '{{if ($3 ~ /unclassified/ || $3 ~ /Eukaryota/) print $2}}' {input} > {output} && date) &> >(tee {log})"

rule extract_reads:
    input:
        ids=rules.headers.output,
        dedup1=rules.deduplicate.output.odup1,
        dedup2=rules.deduplicate.output.odup2
    output:
        ex1=os.path.join(RESULTS_DIR, "euk_bin/extracted/{sample}_R1.fastq"),
        ex2=os.path.join(RESULTS_DIR, "euk_bin/extracted/{sample}_R2.fastq")
    log:
        os.path.join(RESULTS_DIR, "logs/extract_read.euk_bin.{sample}.log")
    conda:
        os.path.join(ENV_DIR, "seqtk.yaml")
    message:
        "Extracting reads from {wildcards.sample}"
    shell:
        "(date && seqtk subseq {input.dedup1} {input.ids} > {output.ex1} && seqtk subseq {input.dedup2} {input.ids} > {output.ex2} && date) &> >(tee {log})"

# TODO
rule concatenate:
    input:
        read1=expand(os.path.join(RESULTS_DIR, "euk_bin/extracted/{sample}_R1.fastq"), sample=SAMPLES),
        read2=expand(os.path.join(RESULTS_DIR, "euk_bin/extracted/{sample}_R2.fastq"), sample=SAMPLES)
    output:
        or1=os.path.join(RESULTS_DIR, "euk_bin/preproc/concat/merged_R1.preprocessed.fastq.gz"),
        or2=os.path.join(RESULTS_DIR, "euk_bin/preproc/concat/merged_R2.preprocessed.fastq.gz")
    log:
        os.path.join(RESULTS_DIR, "logs/concatenate.euk_bin.log")
    message:
        "Concatenating the reads for co-assembly"
    shell:
        "(date && cat {input.read1} > {output.or1} && cat {input.read2} > {output.or2} && date) &> >(tee {log})"

rule deduplicate2:
    input:
        r1=rules.concatenate.output.or1,
        r2=rules.concatenate.output.or2
    output:
        odup1=os.path.join(RESULTS_DIR, "euk_bin/preproc/concat2/merged_R1.preprocessed.fastq.gz"),
        odup2=os.path.join(RESULTS_DIR, "euk_bin/preproc/concat2/merged_R2.preprocessed.fastq.gz")
    log:
        os.path.join(RESULTS_DIR, "logs/deduplicate_round2.euk_bin.log")
    threads:
        config["clumpify"]["threads"]
    conda:
        os.path.join(ENV_DIR, "bbmap.yaml")
    message:
        "Removing duplicate reads for easier downstream assembly"
    shell:
        "(date && clumpify.sh in={input.r1} in2={input.r2} out={output.odup1} out2={output.odup2} dupedist={config[clumpify][dupedist]} dedupe=t optical=t threads={threads} groups={config[clumpify][groups]} -Xmx{config[clumpify][memory]} && date) &> >(tee {log})"


############
# Assembly #
############
rule assembly_sr_megahit:
    input:
        sr1=rules.deduplicate2.output.odup1,
        sr2=rules.deduplicate2.output.odup2
    output:
        os.path.join(RESULTS_DIR, "euk_bin/assembly/ASSEMBLY.fasta")
    log:
        os.path.join(RESULTS_DIR, "logs/megahit.euk_bin.log")
    threads:
        config["megahit"]["threads"]
    conda:
        os.path.join(ENV_DIR, "megahit.yaml")
    message:
        "Assembly: short reads: MEGAHIT"
    shell:
        "(date && megahit -1 {input.sr1} -2 {input.sr2} --kmin-1pass -m 0.9 --k-list 27,37,47,57,67,77,87 --min-contig-len 1000 -t {threads} -o $(dirname {output})/tmp && "
        "cd $(dirname {output}) && "
        "rsync -avP tmp/ . && "
        "ln -sf final.contigs.fa $(basename {output}) && "
        "rm -rf tmp/ && "
        "date) &> >(tee {log})"
            

############
# Analyses #
############
rule eukrep_fasta:
    input:
        rules.assembly_sr_megahit.output
    output:
        os.path.join(RESULTS_DIR, "euk_bin/eukrep/eukrep_contigs.fa")
    log:
        os.path.join(RESULTS_DIR, "logs/eukrep_fasta.euk_bin.log")
    conda:
        os.path.join(ENV_DIR, "eukrep.yaml")
    threads:
        config["threads"]
    shell:
        "(date && EukRep -i {input} -o {output} --min 2000 -m strict && date)"

rule coverm:
    input:
        ref=rules.eukrep_fasta.output,
        read1=expand(os.path.join(RESULTS_DIR, "euk_bin/preproc/concat2/merged_R1.preprocessed.fastq.gz"), sample=SAMPLES),
        read2=expand(os.path.join(RESULTS_DIR, "euk_bin/preproc/concat2/merged_R2.preprocessed.fastq.gz"), sample=SAMPLES)
    output:
        os.path.join(RESULTS_DIR, "euk_bin/coverm/")
    log:
        os.path.join(RESULTS_DIR, "logs/coverm.euk_bin.log")
    threads:
        config["coverm"]["threads"]
    conda:
        os.path.join(ENV_DIR, "coverm.yaml")
    message:
        "Align short reads against eukaryotic contigs: COVERM"
    shell:
        "(date && TMPDIR={RESULTS_DIR} coverm make -o {output} -t {threads} -r {input.ref} -c {input.read1} {input.read2} && date) &> >(tee {log})"

rule concoct_prepare:
    input:
        contigs=rules.eukrep_fasta.output,
	    bam=os.path.join(RESULTS_DIR, "euk_bin/coverm/")
    output:
        coverage=os.path.join(RESULTS_DIR, "euk_bin/concoct/input/concoct_coverage_table.tsv"),
	    contigs_cut=os.path.join(RESULTS_DIR, "euk_bin/concoct/input/contigs_10k.fa")
    log:
        os.path.join(RESULTS_DIR, "logs/concoct_prepare.euk_bin.log")
    conda:
        os.path.join(ENV_DIR, "concoct.yaml")
    message:
        "Metagenomic binning first step: CONCOCT"
    shell:
        """
        cut_up_fasta.py {input.contigs} -c 10000 -o 0 --merge_last -b contigs_10K.bed > {output.contigs_cut}
        concoct_coverage_table.py contigs_10K.bed {input.bam}/*bam > {output.coverage}
        """

rule concoct:
    input:
        contigs_cut=rules.concoct_prepare.output.contigs_cut,
        coverage=rules.concoct_prepare.output.coverage
    output:
        os.path.join(RESULTS_DIR, "euk_bin/concoct/output"),
    log:
        os.path.join(RESULTS_DIR, "logs/concoct.euk_bin.log")
    conda:
        os.path.join(ENV_DIR, "concoct.yaml")
    threads:
        config["concoct"]["threads"]
    message:
        "Metagenomic binning: CONCOCT"
    shell:
        """
        concoct --coverage_file {input.coverage} --composition_file {input.contigs_cut} -t {threads} -b {output}
        """

rule concoct_fasta:
    input:
        clustering=rules.concoct.output,
        contigs=rules.eukrep_fasta.output
    output:
        os.path.join(RESULTS_DIR, "euk_bin/concoct/bins"),
    log:
        os.path.join(RESULTS_DIR, "logs/concoct_fasta.euk_bin.log")
    conda:
        os.path.join(ENV_DIR, "concoct.yaml")
    message:
        "Metagenomic binning: CONCOCT"
    shell:
        """
        merge_cutup_clustering.py {input.clustering}/clustering_gt1000.csv > {input.clustering}/clustering_merged.csv
        extract_fasta_bins.py {input.contigs} {input.clustering}/clustering_merged.csv --output_path {output}        
        """
