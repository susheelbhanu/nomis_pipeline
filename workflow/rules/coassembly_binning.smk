# coassembly_binning workflow
# Runs coassembly of all samples followed by binning of coassembly

localrules: depth

############
# Params
BWA_IDX_EXT = ["amb", "ann", "bwt", "pac", "sa"]            

###########################
# default
rule coassembly_binning: 
    input: 
        DUMMY=os.path.join(RESULTS_DIR, "coassembly/binning/dastool/coassembly.done")
    output:
        touch("status/coassembly_binning.done")


#################
# Preprocessing #
#################
rule coassembly_concat:
    input:
        read1=expand(os.path.join(RESULTS_DIR, "euk_bin/dedup/{sample}_R1.fastq.gz"), sample=SAMPLES),
        read2=expand(os.path.join(RESULTS_DIR, "euk_bin/dedup/{sample}_R2.fastq.gz"), sample=SAMPLES)
    output:
        or1=os.path.join(RESULTS_DIR, "coassembly/concat/merged_R1.preprocessed.fastq.gz"),
        or2=os.path.join(RESULTS_DIR, "coassembly/concat/merged_R2.preprocessed.fastq.gz")
    log:
        os.path.join(RESULTS_DIR, "logs/coassembly_concatenate.log")
    message:
        "Concatenating the reads for co-assembly"
    shell:
        "(date && cat {input.read1} > {output.or1} && cat {input.read2} > {output.or2} && date) &> >(tee {log})"

rule coassembly_dedup:
    input:
        r1=rules.coassembly_concat.output.or1,
        r2=rules.coassembly_concat.output.or2
    output:
        odup1=os.path.join(RESULTS_DIR, "coassembly/dedup/merged_R1.preprocessed.fastq.gz"),
        odup2=os.path.join(RESULTS_DIR, "coassembly/dedup/merged_R2.preprocessed.fastq.gz")
    log:
        os.path.join(RESULTS_DIR, "logs/deduplicate_coassembly_reads.log")
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
rule coassembly_sr_megahit:
    input:
        sr1=rules.coassembly_dedup.output.odup1,
        sr2=rules.coassembly_dedup.output.odup2
    output:
        os.path.join(RESULTS_DIR, "coassembly/assembly/COASSEMBLY.fasta")
    log:
        os.path.join(RESULTS_DIR, "logs/coassembly_megahit.log")
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


############################################
# Coverage and depth estimation
rule coverage:
    input:
        fa=rules.coassembly_sr_megahit.output,
        r1=rules.coassembly_dedup.output.odup1,
        r2=rules.coassembly_dedup.output.odup2
    output:
        os.path.join(RESULTS_DIR, "coassembly/coverage/covstats.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/coassembly.coverage.log")
#    params:
#        tmp_dir=config["coverm"]["tmp_dir"]
    threads:
        config["coverm"]["bigmem_threads"]
    conda:
        os.path.join(ENV_DIR, "coverm.yaml")
    message:
        "calculating coverage for coassembly"
    shell:
        "(date && coverm contig -1 {input.r1} -2 {input.r2} --reference {input.fa} --output-file {output} -t {threads} && date) &> {log}"

rule depth:
    input:
        rules.coverage.output
    output:
        os.path.join(RESULTS_DIR, "coassembly/coverage/coverm_depth.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/coassembly.depth.log")
    message:
        "Getting the depth files required for binning for coassembly"
    shell:
        "(date && tail -n +2 {input} > {output} && date) &> {log}"

# BWA index
rule mapping_index:
    input:
        rules.coassembly_sr_megahit.output
    output:
        expand(os.path.join(RESULTS_DIR, "coassembly/mapping/coassembly.{ext}"), ext=BWA_IDX_EXT)
    log:
        os.path.join(RESULTS_DIR, "logs/mapping.bwa.index.log")
    threads:
        config["bwa"]["threads"]
    params:
        idx_prefix=lambda wildcards, output: os.path.splitext(output[0])[0]
    conda:
        os.path.join(ENV_DIR, "racon.yaml")
    message:
        "Mapping: BWA index for assembly mapping"
    shell:
        "(date && bwa index {input} -p {params.idx_prefix} && date) &> {log}"

# Short reads
rule mapping:
    input:
        r1=rules.coassembly_dedup.output.odup1,
        r2=rules.coassembly_dedup.output.odup2,
        asm=rules.coassembly_sr_megahit.output,
        idx=expand(os.path.join(RESULTS_DIR, "coassembly/mapping/coassembly.{ext}"), ext=BWA_IDX_EXT)
    output:
        os.path.join(RESULTS_DIR, "coassembly/mapping/coassembly.sorted.bam")
    log:
        os.path.join(RESULTS_DIR, "logs/mapping.bwa.mem.coassembly.log")
    threads:
        config["bwa"]["map_threads"]
    params:
        idx_prefix=lambda wildcards, input: os.path.splitext(input.idx[0])[0],
        bam_prefix=lambda wildcards, output: os.path.splitext(output[0])[0],
        chunk_size=config["samtools"]["sort"]["chunk_size"]
    conda:
        os.path.join(ENV_DIR, "racon.yaml")
    message:
        "Mapping short reads to assembly w/ BWA"
    shell:
        "(date && "
        "bwa mem -t {threads} {params.idx_prefix} {input.r1} {input.r2} | "
        "samtools view -@ {threads} -SbT {input.asm} | "
        "samtools sort -@ {threads} -m {params.chunk_size} -T {params.bam_prefix} -o {output} && "
        "date) &> {log}"

rule summarise_depth:
    input:
        rules.mapping.output
    output:
        depth=os.path.join(RESULTS_DIR, "coassembly/coverage/coassembly_depth.txt"),
        paired=os.path.join(RESULTS_DIR, "coassembly/coverage/coassembly_paired.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/summarise_depth.coassembly.log")
    threads:
        config["bwa"]["threads"]
    conda:
        os.path.join(ENV_DIR, "metabat2.yaml")
    message:
        "Getting coverage for coassembly"
    shell:
        "(date && "
        "jgi_summarize_bam_contig_depths --outputDepth {output.depth} --pairedContigs {output.paired} {input} && date) &> {log}"

############################################
# Binning
rule maxbin2:
    input:
        fa=rules.coassembly_sr_megahit.output,
        cov=rules.summarise_depth.output.depth
    output:
        os.path.join(RESULTS_DIR, "coassembly/binning/maxbin2/coassembly.tooshort")
    log:
        os.path.join(RESULTS_DIR, "logs/maxbin2.coassembly.log")
    threads:
        config["maxbin2"]["threads"]
    conda:
        os.path.join(ENV_DIR, "maxbin2.yaml")
    message:
        "Running Maxbin2 on coassembly"
    shell:
        "(date && export PATH=$PATH:{config[maxbin2][perl]} && "
        "run_MaxBin.pl -thread {threads} -contig {input.fa} -out $(dirname {output})/coassembly -abund {input.cov} -min_contig_length {config[maxbin2][min_length]} && date) &> {log}"

rule metabat2:
    input:
        fa=rules.coassembly_sr_megahit.output,
        cov=rules.summarise_depth.output.depth
    output:
        os.path.join(RESULTS_DIR, "coassembly/binning/metabat2/coassembly.tooShort.fa")
    log:
        os.path.join(RESULTS_DIR, "logs/metabat2.coassembly.log")
    threads:
        config["metabat2"]["threads"]
    conda:
        os.path.join(ENV_DIR, "metabat2.yaml")
    message:
        "Running MetaBAT2 on coassembly"
    shell:
        "(date && metabat2 -i {input.fa} -a {input.cov} -o $(dirname {output})/coassembly -t {threads} -m {config[metabat2][min_length]} -v --unbinned --cvExt && date) &> {log}"

############################################
# DASTool
rule scaffold:
    input:
        max=rules.maxbin2.output,
        met=rules.metabat2.output
    output:
        maxscaf=os.path.join(RESULTS_DIR, "coassembly/binning/maxbin2/coassembly_maxbin2.scaffolds2bin.tsv"),
        metscaf=os.path.join(RESULTS_DIR, "coassembly/binning/metabat2/coassembly_metabat2.scaffolds2bin.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/scaffold.coassembly.log")
    conda:
        os.path.join(ENV_DIR, "dastool.yaml")
    message:
        "Getting DASTool:scaffold list for coassembly"
    shell:
        "(date && scripts/Fasta_to_Scaffolds2Bin.sh -i $(dirname {input.max}) -e fa > {output.maxscaf} && "
        "scripts/Fasta_to_Scaffolds2Bin.sh -i $(dirname {input.met}) -e fa > {output.metscaf} && date) &> {log}"

rule dastool:
    input:
        max=rules.scaffold.output.maxscaf,
        met=rules.scaffold.output.metscaf,
        fa=rules.coassembly_sr_megahit.output
    output:
        DIR=directory(os.path.join(RESULTS_DIR, "coassembly/binning/dastool/")),
        DUMMY=os.path.join(RESULTS_DIR, "coassembly/binning/dastool/coassembly.done")
    log:
        os.path.join(RESULTS_DIR, "logs/dastool.coassembly.log")
    threads:
        config["dastool"]["threads"]
    conda:
        os.path.join(ENV_DIR, "dastool.yaml")
    params:
        db=config["dastool"]["database"]
    message:
        "Running DASTool binning on coassembly"
    shell:
        "(date && export PATH=$PATH:{config[dastool][path]} && "
        "export PATH=$PATH:{config[dastool][src]} && "
        "DAS_Tool -i {input.max},{input.met} -c {input.fa} -o $(dirname {output.DIR}) --score_threshold {config[dastool][score]} --search_engine diamond -l maxbin2,metabat2 --write_bins 1 --write_bin_evals 1 --threads {threads} --db_directory {params.db} --create_plots 1 && "
        "touch {output.DUMMY} && date) &> {log}"
