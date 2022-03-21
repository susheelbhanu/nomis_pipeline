# imp workflow
# Sets up folders and IMP launch scripts per sample

localrules: folders, imp_files

###########################
# default

rule imp: 
    input: 
        expand(os.path.join(DATA_DIR, "{sample}/{sample}_Pool1_R{read}.fastq.gz"), read=["1", "2"], sample=SAMPLES),
        expand(os.path.join(DATA_DIR, "{sample}/{sample}_config.yaml"), sample=SAMPLES),
        expand(os.path.join(DATA_DIR, "{sample}/launchIMP.sh"), sample=SAMPLES),
        expand(os.path.join(DATA_DIR, "{sample}/{sample}.runIMP.sh"), sample=SAMPLES)
    output:
        touch("status/imp.done")


################################
# rules for files and analyses #
################################
rule folders:
    input:
        in1=os.path.join(FASTQ_DIR, "20200624.A-{sample}_Pool1_R1.fastq.gz"),
        in2=os.path.join(FASTQ_DIR, "20200624.A-{sample}_Pool1_R2.fastq.gz")
    output:
        fout1=os.path.join(DATA_DIR, "{sample}/{sample}_Pool1_R1.fastq.gz"),
        fout2=os.path.join(DATA_DIR, "{sample}/{sample}_Pool1_R2.fastq.gz")
    shell:
        "(date && ln -vs {input.in1} {output.fout1} && "
        "ln -vs {input.in2} {output.fout2} && date) &> >(tee {log})"

rule imp_files:
    input:
        config=os.path.join(NOTES_DIR, "IMP_config.yaml"),
        runfile=os.path.join(NOTES_DIR, "runIMP.sh"),
        launcher=os.path.join(NOTES_DIR, "launchIMP.sh")
    output:
        tout1=os.path.join(DATA_DIR, "{sample}/{sample}_config.yaml"),
        tout2=os.path.join(DATA_DIR, "{sample}/launchIMP.sh"),
        tout3=os.path.join(DATA_DIR, "{sample}/{sample}.runIMP.sh")
    shell:
        "(date && cp -v {input.config} {output.tout1} && "
        "cp -v {input.launcher} {output.tout2} && "
        "cp -v {input.runfile} {output.tout3} && "
        "sed -i 's/\"\$sample\"/{wildcards.sample}/g' {output.tout1} && "
        "sed -i 's/\"\$sample\"/{wildcards.sample}/g' {output.tout2} && "
        "sed -i 's/\"\$sample\"/{wildcards.sample}/g' {output.tout3} && date) &> >(tee {log})"