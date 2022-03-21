# viruses workflow
# Runs VIBRANT, vCONTACT2, Kaiju with viral database and CheckV

###########################
# default

rule viruses: 
    input: 
        expand(os.path.join(RESULTS_DIR, "vibrant_output/VIBRANT_{sample}/{sample}.prodigal.faa"), sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "vibrant_output/VIBRANT_{sample}/VIBRANT_phages_{sample}/{sample}.phages_combined.simple.faa"), sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "vibrant_output/VIBRANT_{sample}/{sample}.g2g.csv"), sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "vcontact2_output/{sample}/viral_cluster_overview.csv"), sample=SAMPLES),
#        expand(os.path.join(RESULTS_DIR, "vcontact2_output/{sample}_{file}.csv"), sample=SAMPLES, file=["viral_cluster_overview", "genome_by_genome_overview", "merged_df", "vConTACT_contigs"]),
        expand(os.path.join(RESULTS_DIR, "kraken2_output/{sample}_kraken.report.txt"), sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "kaiju_output/{sample}_kaiju.out"), sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "kaiju_output/kaiju_summary.out")),
        expand(os.path.join(RESULTS_DIR, "checkv/{sample}/quality_summary.tsv"), sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "antismash/{sample}/knownclusterblastoutput.txt"), sample=SAMPLES)
    output:
        touch("status/viruses.done")


################################
# rules for files and analyses #
################################
rule vibrant:
    input:
        os.path.join(DATA_DIR, "{sample}/run1/Assembly/mg.assembly.merged.fa")
    output:    
#        viout1=os.path.join(DATA_DIR, "vibrant_output/VIBRANT_{sample}/VIBRANT_log_annotation_{sample}.log")
#        viout2=directory(os.path.join(DATA_DIR, "vibrant_output/"))
        viout1=os.path.join(RESULTS_DIR, "vibrant_output/VIBRANT_{sample}/{sample}.prodigal.faa"),
        viout2=os.path.join(RESULTS_DIR, "vibrant_output/VIBRANT_{sample}/VIBRANT_phages_{sample}/{sample}.phages_combined.faa")
    threads:
        config['threads']
    conda:
        os.path.join(ENV_DIR, "vibrant.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/vibrant.{sample}.log")
    message:
        "Running vibrant for {wildcards.sample}"
    shell:
        "(date && python3 ./vibrant/VIBRANT/VIBRANT_run.py -t {threads} -i {input} -folder $(dirname $(dirname {output.viout1})) && date) &> >(tee {log})"

rule convert_files:
    input:
        os.path.join(RESULTS_DIR, "vibrant_output/VIBRANT_{sample}/VIBRANT_phages_{sample}/{sample}.phages_combined.faa")
    output:
        tout1=os.path.join(RESULTS_DIR, "vibrant_output/VIBRANT_{sample}/VIBRANT_phages_{sample}/{sample}.phages_combined.simple.faa"),
        tout2=os.path.join(RESULTS_DIR, "vibrant_output/VIBRANT_{sample}/{sample}.g2g.csv")
    threads:
        config['convert_files']['threads']
#    conda:
#        os.path.join(ENV_DIR, "vcontact2.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/convert.{sample}.log")
    message:
        "Converting files for {wildcards.sample}"
    shell:
        "(date && python3 {config[convert_files][simplify]} {input} && "
        "export PATH='/scratch/users/sbusi/tools/miniconda3/envs/vcontact2/bin:$PATH' && "
        "vcontact2_gene2genome -p {output.tout1} -o {output.tout2} -s '{config[convert_files][type]}') &> >(tee {log})"

rule vcontact2:
    input:
        v1=rules.convert_files.output.tout1,
        v2=rules.convert_files.output.tout2
    output:
        cout1=os.path.join(RESULTS_DIR, "vcontact2_output/{sample}/viral_cluster_overview.csv"),
        cout2=os.path.join(RESULTS_DIR, "vcontact2_output/{sample}/genome_by_genome_overview.csv"),
        cout3=os.path.join(RESULTS_DIR, "vcontact2_output/{sample}/merged_df.csv"),
        cout4=os.path.join(RESULTS_DIR, "vcontact2_output/{sample}/vConTACT_contigs.csv"),
        cout5=directory(os.path.join(RESULTS_DIR, "vcontact2_output/{sample}"))
    threads:
        config['convert_files']['threads']
#    conda:
#        os.path.join(ENV_DIR, "vcontact2.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/vcontact2.{sample}.log")
    message:
        "Running vcontact2 on {wildcards.sample}"
    shell:
        "(date && export PATH=$PATH:'/scratch/users/sbusi/tools/miniconda3/envs/vcontact2/bin' && "
        "/scratch/users/sbusi/tools/miniconda3/envs/vcontact2/bin/vcontact2 --force-overwrite --raw-proteins {input.v1} --rel-mode 'Diamond' --proteins-fp {input.v2} --db 'ProkaryoticViralRefSeq94-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin /home/users/sbusi/apps/vcontact2/cluster_one-1.0.jar --output-dir {output.cout5} && date) &> >(tee {log})"

# rule copy:
#     input: 
#         cop1=os.path.join(RESULTS_DIR, "vcontact2_output/{sample}/viral_cluster_overview.csv"),
#         cop2=os.path.join(RESULTS_DIR, "vcontact2_output/{sample}/genome_by_genome_overview.csv"),
#         cop3=os.path.join(RESULTS_DIR, "vcontact2_output/{sample}/merged_df.csv"),
#         cop4=os.path.join(RESULTS_DIR, "vcontact2_output/{sample}/vConTACT_contigs.csv")
#     output:
#         vout1=os.path.join(RESULTS_DIR, "vcontact2_output/{sample}_viral_cluster_overview.csv"),
#         vout2=os.path.join(RESULTS_DIR, "vcontact2_output/{sample}_genome_by_genome_overview.csv"),
#         vout3=os.path.join(RESULTS_DIR, "vcontact2_output/{sample}_merged_df.csv"),
#         vout4=os.path.join(RESULTS_DIR, "vcontact2_output/{sample}_vConTACT_contigs.csv")
#     log:
#         os.path.join(RESULTS_DIR, "logs/copy.{sample}.log")
#     message:
#         "Copying files from {wildcards.sample}"
#     shell:
#         "(date && ln -vs {input.cop1} {output.vout1} && ln -vs {input.cop2} {output.vout2} && ln -vs {input.cop3} {output.vout3} && ln -vs {input.cop4} {output.vout4} && date) &> >(tee {log})"

##################
# Viral taxonomy #
##################
rule kraken2:
    input:
        os.path.join(DATA_DIR, "{sample}/run1/Assembly/mg.assembly.merged.fa")
    output:
        report=os.path.join(RESULTS_DIR, "kraken2_output/{sample}_kraken.report.txt"),
        summary=os.path.join(RESULTS_DIR, "kraken2_output/{sample}_kraken.out")
    conda:
        os.path.join(ENV_DIR, "kraken2.yaml")
    threads:
        config['kraken2']['threads']
    log:
        os.path.join(RESULTS_DIR, "logs/kraken2.{sample}.log")
    message:
        "Running kraken2 on {wildcards.sample}"
    shell:
        "(date && kraken2 --threads {threads} --db {config[kraken2][db]} --confidence 0.75 {input} --output {output.summary} --report {output.report} && date) &> >(tee {log})"

rule kaiju:
    input:
        fasta=os.path.join(DATA_DIR, "{sample}/run1/Assembly/mg.assembly.merged.fa"),
#        nodes=os.path.join(DATA_DIR, "kaiju/db/viruses/nodes.dmp"),
#        fi=os.path.join(DATA_DIR, "kaiju/db/viruses/kaiju_db_viruses.fmi")
    output:
        os.path.join(RESULTS_DIR, "kaiju_output/{sample}_kaiju.out")
    params:
        nodes="nodes.dmp",
        fmi="kaiju_db_viruses.fmi"
    conda:
        os.path.join(ENV_DIR, "kaiju.yaml")
    threads:
        config['kaiju']['threads']
    log:
        os.path.join(RESULTS_DIR, "logs/kaiju.{sample}.log")
    message:
        "Tax. classification w/ Kaiju viral-db for {wildcards.sample}"
    shell:
        "(date && kaiju -z {threads} -t {config[kaiju][db]}/{params.nodes} -f {config[kaiju][db]}/{params.fmi} -i {input.fasta} -o {output} && date) &> >(tee {log})"

rule kaiju_summary:
    input:
        files=expand(os.path.join(RESULTS_DIR, "kaiju_output/{sample}_kaiju.out"),  sample=SAMPLES),
        # nodes=os.path.join(DATA_DIR, "kaiju/db/viruses/nodes.dmp"),
        # names=os.path.join(DATA_DIR, "kaiju/db/viruses/names.dmp")
    output:
        os.path.join(RESULTS_DIR, "kaiju_output/kaiju_summary.out")
    params:
        nodes="nodes.dmp",
        names="names.dmp"
    conda:
        os.path.join(ENV_DIR, "kaiju.yaml")
    threads: 1
    log:
        os.path.join(RESULTS_DIR, "logs/kaiju-summary.log")
    message:
        "Summarizing Kaiju viral-db output"
    shell:
        "(date && kaiju2table -e -t {config[kaiju][db]}/{params.nodes} -n {config[kaiju][db]}/{params.names} -r {config[kaiju][rank]} -o {output} {input.files} && date) &> >(tee {log})"

rule checkv:
    input:
        os.path.join(DATA_DIR, "{sample}/run1/Assembly/mg.assembly.merged.fa")
    output:
        os.path.join(RESULTS_DIR, "checkv/{sample}/quality_summary.tsv")
    conda:
        os.path.join(ENV_DIR, "checkv.yaml")
    threads:
        config['checkv']['threads']
    log:
        os.path.join(RESULTS_DIR, "logs/checkv.{sample}.log")
    message:
        "Running CheckV for {wildcards.sample}"
    shell:
        "(date && checkv end_to_end -d {config[checkv][db]} {input} $(dirname {output}) -t {threads} && date) &> >(tee {log})"

######################
# rule for antiSMASH #
######################
rule antismash:
    input:
        FA=os.path.join(DATA_DIR, "{sample}/run1/Assembly/mg.assembly.merged.fa"),
        GFF=os.path.join(DATA_DIR, "{sample}/run1/Analysis/annotation/annotation_CDS_RNA_hmms.gff")
    output:
        os.path.join(RESULTS_DIR, "antismash/{sample}/knownclusterblastoutput.txt")
    conda:
        os.path.join(ENV_DIR, "antismash.yaml")
    threads:
        config['antismash']['threads']
    log:
        os.path.join(RESULTS_DIR, "logs/antismash.{sample}.log")
    message:
        "Running antiSMASH for {wildcards.sample}"
    shell:
        "(date && antismash --cpus {threads} --genefinding-tool none --genefinding-gff3 {input.GFF} --fullhmmer --pfam2go --asf --cb-knownclusters --clusterhmmer --cf-create-clusters {input.FA} --output-dir $(dirname {output}) && date) &> >(tee {log})"