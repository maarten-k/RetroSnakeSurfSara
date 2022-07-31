configfile: "config.yaml"


SAMPLES = ["samplecram"]
outPath = config["outPath"]
bamPath = config["bamPath"]
cramPath = config["cramPath"]


containerized: "/home/maartenk/temp/RetroSnakeSurfSara/sif/sif/retroseq.sif"


localrules:
    install_modified_retroseq,


rule all:
    input:
        expand(outPath + "filter/{sample}.bed", sample=SAMPLES),
        expand(outPath + "confirmed/{sample}.retroseqHitsConfirmed.bed", sample=SAMPLES),


rule install_modified_retroseq:
    output:
        "resources/RetroSeq/bin/retroseq.pl",
    threads: 1
    shell:
        """
        mkdir -p resources/RetroSeq/
        cd resources/RetroSeq/
        wget https://github.com/maarten-k/RetroSeq/archive/refs/heads/master.zip
        unzip master.zip
        rm master.zip
        mv RetroSeq-master/* .
        rmdir RetroSeq-master
        """


rule retroseqDiscover:
    input:
        cram=cramPath + "{sample}.cram",
        index=cramPath + "{sample}.cram.crai",
        retroseq="resources/RetroSeq/bin/retroseq.pl",
    output:
        outPath + "discover/{sample}.bed",
    threads: 8
    log:
        "logs/discover/{sample}.log",
    params:
        identity=80,
    benchmark:
        #repeat("benchmarks/{sample}.retroseqDiscover.benchmark.txt",3)
        "benchmarks/{sample}.retroseqDiscover.benchmark.txt"
    conda:
        "envs/retroseq.yaml"
    shell:
        "perl {input.retroseq} -discover -bam {input.cram} -output {output} -eref {config[HERVK_eref]} -id {params.identity}"


rule retroseqCall:
    input:
        bam=cramPath + "{sample}.cram",
        bai=cramPath + "{sample}.cram.crai",
        discover=outPath + "discover/{sample}.bed",
        retroseq="resources/RetroSeq/bin/retroseq.pl",
    output:
        outPath + "call/{sample}.vcf",
    threads: 8
    benchmark:
        "benchmarks/{sample}.retroCall.benchmark.txt"
    conda:
        "envs/retroseq.yaml"
    log:
        "logs/call/{sample}.log",
    shell:
        "{input.retroseq} -call -bam {input.bam} -input {input.discover} -ref {config[refHg19]} -output {output}"


rule filterCalls:
    input:
        outPath + "call/{sample}.vcf",
    output:
        outPath + "filter/{sample}.pos",
        outPath + "filter/{sample}.bed",
    log:
        "logs/filter/{sample}.log",
    conda:
        "envs/verification.yaml"
    benchmark:
        "benchmarks/{sample}.filterCalls.benchmark.txt"
    shell:
        "python {config[pythonScripts]}/filterHighQualRetroseqForDownstream.py  {input} {output}"


rule verify:
    input:
        outPath + "filter/{sample}.pos",
        cramPath + "{sample}.cram",
        cramPath + "{sample}.cram.crai",
    output:
        outPath + "confirmed/{sample}.retroseqHitsConfirmed.bed",
    benchmark:
        "benchmarks/{sample}.verify.tsv"
    params:
        verificationLevel="low",
    conda:
        "envs/verification.yaml"
    log:
        "logs/call/{sample}.log",
    shell:
        """
        python {config[pythonScripts]}/assembleAndRepeatMasker.py {input[0]} {config[cramPath]}{wildcards.sample}.cram {config[outPath]} $CONDA_PREFIX/bin/ {config[pythonScripts]} {config[element]} {params.verificationLevel} {output} 
        """
