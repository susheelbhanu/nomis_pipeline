# functional analyses workflow
# Running functional analyses tools: METABOLIC + MAGICCAVE + MetabolisHMM + MANTIS + FUNCS

localrules: metadata, install_magiccave

###########################
# default

rule functions:
    input:
        os.path.join(RESULTS_DIR, "metabolic_output"),
        expand(os.path.join(RESULTS_DIR, "funcs/gtdbtk/{sample}_gtdbtk.txt"), sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "funcs/kegg/{sample}_kegg.txt"), sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "funcs/kegg/{sample}_cov_len.txt"), sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "funcs/scaffolds/{sample}_scaffolds2bin.txt"), sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "funcs/merged_files/merged_{filetype}.txt"), filetype=["gtdbtk", "kegg", "scaffolds2bin", "kegg_cov_len"]),
        os.path.join(RESULTS_DIR, "funcs/merged_files/ROCK_kegg_taxa_mag.csv"),
        os.path.join(RESULTS_DIR, "metabolishmm/summarize_metabolism_out"),
        os.path.join(RESULTS_DIR, "magiccave/lithogenie_output")
    output:
        touch("status/functions.done")


#######################
# Functional analyses #
#######################
# METABOLIC potential of bins
# NOTE: some adjustments to perl environments made based on here: https://github.com/AnantharamanLab/METABOLIC/issues/27
rule prep_metabolic:
    input:
        os.path.join(RESULTS_DIR, "bins/bin_collection.done")
    output:
        sample=temp(os.path.join(RESULTS_DIR, "data/metabolic_samples.txt")),
        tmpfile=temp(os.path.join(RESULTS_DIR, "data/reads.txt")),
        reads=os.path.join(RESULTS_DIR, "data/metabolic_reads.txt")
    message:
        "Creating read mapping file for METABOLIC"
    shell:
        "(date && "
        "while read -r line; do ls $(dirname {input})/*.fa | grep -o \"$line\" ; done < {output.sample} > {output.tmpfile} && "
        "sed 's@^@/work/projects/nomis/metaG_JULY_2020/IMP3/@g' {output.tmpfile} | "
        "sed 's@$@/run1/Preprocessing/mg.r1.preprocessed.fq@g' | "
        "awk -F, '{{print $0=$1\",\"$1}}' | awk 'BEGIN{{FS=OFS=\",\"}} {{gsub(\"r1\", \"r2\", $2)}} 1' | "
        "sed $'1 i\\\\\\n# Read pairs:' {output.reads}"     # using forward-slashes to get `\\\n`

rule metabolic:
    input:
        fa=os.path.join(RESULTS_DIR, "bins/bin_collection.done"),
        reads=rules.prep_metabolic.output
    output:
        directory(os.path.join(RESULTS_DIR, "metabolic_output"))
    log:
        os.path.join(RESULTS_DIR, "logs/metabolic.log")
    conda:
        os.path.join(ENV_DIR, "metabolic.yaml")
    params:
        gtdbtk=config["metabolic"]["db"],
        metabolic=config["metabolic"]["directory"]
    threads:
        config["metabolic"]["threads"]
    message:
        "Running metabolic for all MAGs"
    shell:
        "(date && "
        "export GTDBTK_DATA_PATH={params.gtdbtk} && "
        "export PERL5LIB && export PERL_LOCAL_LIB_ROOT && export PERL_MB_OPT && export PERL_MM_OPT && "
        """env PERL5LIB="" PERL_LOCAL_LIB_ROOT="" PERL_MM_OPT="" PERL_MB_OPT="" cpanm Array::Split && """
        "perl {params.metabolic}/METABOLIC-C.pl -t {threads} -in-gn $(dirname {input.fa}) -r {input.reads} -o {output} && "
        "date) &> >(tee {log})"


###########################################
# FUNCS analyses: rules for file creation #
###########################################
rule func_bins:
    input:
        os.path.join(DATA_DIR, "{sample}/run1/Stats/all_bin_stats.tsv")
    output:
        os.path.join(RESULTS_DIR, "funcs/gtdbtk/{sample}_gtdbtk.txt")
    log:
        os.path.join(DATA_DIR, "logs/bins.{sample}.log")
    message:
        "Getting the bins-taxa file for {wildcards.sample}"
    run:
        bin=pd.read_csv(input[0], sep="\t")
        bin_edited=bin[['selected_by_DASTool', 'classification']]       # selecting columns
#        bins=bin_edited[~bin_edited.selected_by_DASTool.str.contains("NaN", na=False)]  # dropping rows where bins not selected by DASTool
        bins=bin_edited.dropna(subset=['selected_by_DASTool'])
        bins.rename(columns = {'selected_by_DASTool': 'Bin', 'classification': 'Taxa'}, inplace=True)
        bins.to_csv(output[0], sep="\t", index=False)

rule kegg_contigs:
    input:
        os.path.join(DATA_DIR, "{sample}/run1/Analysis/annotation/mg.KEGG.counts.tsv")
    output:
        os.path.join(RESULTS_DIR, "funcs/kegg/{sample}_kegg.txt")
    log:
        os.path.join(DATA_DIR, "logs/kegg.{sample}")
    message:
        "Getting the kegg-contigs file for {wildcards.sample}"
    run:
        kegg=pd.read_csv(input[0], sep="\t", skiprows=1)
        kegg_edited=kegg[['Geneid', 'Chr']]
        kegg_edited.rename(columns = {'Geneid': 'KEGG', 'Chr': 'Contig'}, inplace=True)
        kegg_contigs=(kegg_edited.assign(Contig = kegg_edited['Contig'].str.split(';')).explode('Contig').reset_index(drop=True))
        kegg_contigs=kegg_contigs.reindex(['Contig','KEGG'], axis=1)
        kegg_contigs.to_csv(output[0], sep="\t", index=False)

rule kegg_coverage_length:
    input:
        kegg=rules.kegg_contigs.output,
        cov=os.path.join(DATA_DIR, "{sample}/run1/Stats/mg/mg.assembly.contig_depth.txt"),
        length=os.path.join(DATA_DIR, "{sample}/run1/Stats/mg/mg.assembly.length.txt")
    output:
        os.path.join(RESULTS_DIR, "funcs/kegg/{sample}_cov_len.txt")
    log:
        os.path.join(DATA_DIR, "logs/kegg_cov_len.{sample}")
    message:
        "Getting the kegg-contigs file for {wildcards.sample}"
    run:
        kegg=pd.read_csv(input[0], sep="\t", header=0)
        cov=pd.read_csv(input[1], sep="\t", header=None)
        cov.rename(columns={0: 'Contig', 1: 'Coverage'}, inplace=True)
        
        length=pd.read_csv(input[2], sep="\t", header=None)
        length.rename(columns={0: 'Contig', 1: 'Length'}, inplace=True)

        tmp=pd.merge(kegg, cov, on='Contig')
        all_merged=pd.merge(tmp, length, on='Contig')
        all_merged.to_csv(output[0], sep="\t", index=False)

rule scaffolds:
    input:
        os.path.join(DATA_DIR, "{sample}/run1/Binning/selected_DASTool_scaffolds2bin.txt")
    output:
        os.path.join(RESULTS_DIR, "funcs/scaffolds/{sample}_scaffolds2bin.txt")
    log: 
        os.path.join(DATA_DIR, "logs/scaffolds.{sample}.log")
    message:
        "Getting the scaffolds2bin file for {wildcards.sample}"
    run:
        scaffold=pd.read_csv(input[0], sep="\t", header=None)
        scaffold.columns=['Contig', 'Bin']
        scaffold.to_csv(output[0], sep="\t", index=False)

# FUNCS: Rules for merging 
rule merge_bins:
    input:
        expand(os.path.join(RESULTS_DIR, "funcs/gtdbtk/{sample}_gtdbtk.txt"), sample=SAMPLES)
    output:
        os.path.join(RESULTS_DIR, "funcs/merged_files/merged_gtdbtk.txt")
    message:
        "Merging all the bin-taxa files into one file"
    run:
        opened=[]
        for ifile in input:
          df=pd.read_csv(ifile, index_col=None, sep="\t", header=0)
          df['Sample']=re.sub("_gtdbtk.txt", "", os.path.basename(ifile))
          df = df.reindex(['Sample','Bin','Taxa'], axis=1)
          opened.append(df)

        frame=pd.concat(opened, axis=0, ignore_index=True)
        frame.to_csv(output[0], sep="\t", index=False)

rule merge_keggs:
    input:
        expand(os.path.join(RESULTS_DIR, "funcs/kegg/{sample}_kegg.txt"), sample=SAMPLES)
    output:
        os.path.join(RESULTS_DIR, "funcs/merged_files/merged_kegg.txt")
    message:
        "Merging all the KEGG-contigs file into one files"
    run:
        opened=[]
        for ifile in input:
          df=pd.read_csv(ifile, index_col=None, sep="\t", header=0)
#          df['Sample']=re.sub("_gtdbtk.txt", "", os.path.basename(ifile))
#          df = df.reindex(['Sample','Contig','KEGG'], axis=1)
          opened.append(df)

        frame=pd.concat(opened, axis=0, ignore_index=True)
        frame.to_csv(output[0], sep="\t", index=False) 

rule merge_keggs_cov_len:
    input:
        expand(os.path.join(RESULTS_DIR, "funcs/kegg/{sample}_cov_len.txt"), sample=SAMPLES)
    output:
        os.path.join(RESULTS_DIR, "funcs/merged_files/merged_kegg_cov_len.txt")
    message:
        "Merging all the files with KEGG-contigs plus coverage and length included into one file"
    run:
        opened=[]
        for ifile in input:
          df=pd.read_csv(ifile, index_col=None, sep="\t", header=0)
          opened.append(df)

        frame=pd.concat(opened, axis=0, ignore_index=True)
        frame.to_csv(output[0], sep="\t", index=False)

rule merge_scaffolds:
    input:
        expand(os.path.join(RESULTS_DIR, "funcs/scaffolds/{sample}_scaffolds2bin.txt"), sample=SAMPLES)
    output:
        os.path.join(RESULTS_DIR, "funcs/merged_files/merged_scaffolds2bin.txt")
    message:
        "Merging all the scaffolds2bin files into one file"
    run:
        opened=[]
        for ifile in input:
          df=pd.read_csv(ifile, index_col=None, sep="\t", header=0)
#          df['Sample']=re.sub("_gtdbtk.txt", "", os.path.basename(ifile))
#          df = df.reindex(['Sample','Bin','Taxa'], axis=1)
          opened.append(df)

        frame=pd.concat(opened, axis=0, ignore_index=True)
        frame.to_csv(output[0], sep="\t", index=False)         

rule merge_all:
    input:
        coverage=rules.merge_keggs_cov_len.output,
        scaffolds=rules.merge_scaffolds.output,
        gtdbtk=rules.merge_bins.output
    output:
        merged=os.path.join(RESULTS_DIR, "funcs/merged_files/ROCK_kegg_taxa_mag.csv")
    conda:
        os.path.join(ENV_DIR, "r-conda.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/merge_all.log")
    message:
        "Merging all kegg, coverage, length, bin, sample, and taxonomy for FUNCs analyses"
    script:
        "merge_funcs.R"


######################################
# rules for metabolisHMM + MagicCave #
######################################
rule metadata:
    input:
        stat=os.path.join(RESULTS_DIR, "gtdbtk_output")
    params:
        arc=os.path.join(RESULTS_DIR, "gtdbtk_output/gtdbtk.ar122.summary.tsv"),
        bac=os.path.join(RESULTS_DIR, "gtdbtk_output/gtdbtk.bac120.summary.tsv")
    output:
        os.path.join(RESULTS_DIR, "data/metadata")
    shell:
        "cat {params} | awk '{{print $1\",\"$2}}' | sed '@^user@d' | sed 's@.fasta.contigs@@g' | sed 's@.fasta_sub.contigs@@g' | awk '!visited[$0]++' > {output}"

rule metabolishmm:
    input:
        bins=os.path.join(RESULTS_DIR, "bins/bin_collection.done"),
        meta=rules.metadata.output
    output:
        sum=directory(os.path.join(RESULTS_DIR, "metabolishmm/summarize_metabolism_out")),
        indiv=directory(os.path.join(RESULTS_DIR, "metabolishmm/individual_metabolism_out"))
    conda:
        os.path.join(ENV_DIR, "metabolishmm.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/metabolishmm.log")
    message:
        "Running MetabolisHMM on all MAGs"
    shell:
        "(date && summarize-metabolism --input $(dirname {input.bins}) --output {output.sum} --metadata {input.meta} --summary {output.sum}/results/summarize_metabolism.csv --heatmap {output.sum}/results/summarize_metabolism.pdf --aggregate ON --plotting ON && summarize-metabolism --input $(dirname {input.bins}) --output {output.indiv} --metadata {input.meta} --summary {output.indiv}/results/individual_metabolism.csv --heatmap {output.indiv}/results/individual_metabolism.pdf --plotting ON && date) &> >(tee {log})"


# MagicLamp Initial Setup 
rule install_magiccave:
    output:
        done=os.path.join(RESULTS_DIR, "magiccave/magiccave.installed")
    log:
        out=os.path.join(RESULTS_DIR, "logs/setup.magiccave.log")
    params:
        script=os.path.join(SRC_DIR, "install_magiccave.sh"), 
        path=os.path.join(SUBMODULES, "magiccave")
    conda:
        os.path.join(ENV_DIR, "magiccave.yaml")
    message:
        "Setup: install MagicCave"
    shell:
        "script=$(realpath {params.script}) && cd {params.path} && ${{script}}"

rule magiccave:
    input:
        bins=os.path.join(RESULTS_DIR, "bins/bin_collection.done"),
        installed=os.path.join(RESULTS_DIR, "magiccave/magiccave.installed")
    output:
        directory(os.path.join(RESULTS_DIR, "magiccave/lithogenie_output"))
    params:
        path=os.path.join(SUBMODULES, "magiccave")
    conda:
        os.path.join(ENV_DIR, "magiccave.yaml")
    threads:
        config['magiccave']['threads']
    log:
        os.path.join(RESULTS_DIR, "logs/magiccave.log")
    message:
        "Running MagicLamps from MagicCave on all MAGs"
    shell:
        "(date && export PATH=$PATH:{params.path} && "
        "MagicLamp.py LithoGenie -bin_dir $(dirname {input.bins}) -bin_ext fa -out {output} -t {threads} --norm && date) &> >(tee {log})" 
