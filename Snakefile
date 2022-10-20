configfile: "config.yaml"


SAMPLES = [config["sample"]
]
outPath = config["outPath"]
cramPath = config["cramPath"]


containerized: config["containerized_path"]


localrules:
    install_modified_retroseq,


rule all:
    input:
        expand(outPath + "filter/{sample}.bed", sample=SAMPLES),
        expand(outPath + "marked/{sample}.novelHitsFV.bed", sample=SAMPLES),
        expand(outPath + "confirmed/{sample}.retroseqHitsConfirmed.bed", sample=SAMPLES),


rule install_modified_retroseq:
    output:
        "resources/RetroSeq/bin/retroseq.pl",
    threads: 1
    container:
        None
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
    threads: 1
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
        "perl {input.retroseq} -discover -bam {input.cram} -output {output} -eref {config[HERVK_eref]} -id {params.identity} > {log}"


rule retroseqCall:
    input:
        bam=cramPath + "{sample}.cram",
        bai=cramPath + "{sample}.cram.crai",
        discover=outPath + "discover/{sample}.bed",
        retroseq="resources/RetroSeq/bin/retroseq.pl",
    output:
        temp(outPath + "call/{sample}_{chr}.vcf"),
    threads: 1
    benchmark:
        "benchmarks/{sample}_{chr}.retroCall.benchmark.txt"
    conda:
        "envs/retroseq.yaml"
    log:
        "logs/call/{sample}_{chr}.log",
    shell:
        "{input.retroseq} -call -region {wildcards.chr} -bam {input.bam} -input {input.discover} -ref {config[refHg19]} -output {output} > {log}"


rule MergeCalls:
    input:
        vcf=expand(
            outPath + "call/{{sample}}_{chr}.vcf",
            chr=[
                "chr1",
                "chr2",
                "chr3",
                "chr4",
                "chr5",
                "chr6",
                "chr7",
                "chr8",
                "chr9",
                "chr10",
                "chr11",
                "chr12",
                "chr13",
                "chr14",
                "chr15",
                "chr16",
                "chr17",
                "chr18",
                "chr19",
                "chr20",
                "chr21",
                "chr22",
                "chrY",
                "chrX",
                "chrM",
            ],
        ),
    output:
        outPath + "call/{sample}.vcf",
    conda:
        "envs/retroseq.yaml"
    threads: 1
    shell:
        """
        grep "^#" {input.vcf[0]} >{output}
        grep -hv "^#" {input.vcf}>>{output}
        """


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


rule markAll:
    input:
        outPath + "filter/{sample}.bed",
        outPath + "confirmed/{sample}.retroseqHitsConfirmed.bed",
    output:
        outPath + "marked/{sample}.knownHitsF.bed",
        outPath + "marked/{sample}.novelHitsF.bed",
        outPath + "marked/{sample}.knownHitsFV.bed",
        outPath + "marked/{sample}.novelHitsFV.bed",
    conda:
        "envs/retroseq.yaml"
    shell:
        """
        bedtools window -w 500 -c -a {config[knownNR]} -b  {input[0]}  > {output[0]}
        bedtools sort -i  {input[0]} >  {config[outPath]}marked/{wildcards.sample}.sorted.bed
        bedtools window -w 500 -v -a {config[outPath]}marked/{wildcards.sample}.sorted.bed -b {config[knownNR]} > {output[1]} 
        rm {config[outPath]}marked/{wildcards.sample}.sorted.bed
        bedtools window -w 500 -c -a {config[knownNR]} -b  {input[1]} > {output[2]}
        bedtools sort -i  {input[1]} >  {config[outPath]}marked/{wildcards.sample}.sorted.bed
        bedtools window -w 500 -v -a {config[outPath]}marked/{wildcards.sample}.sorted.bed -b {config[knownNR]} > {output[3]} 
        rm {config[outPath]}marked/{wildcards.sample}.sorted.bed
        """
