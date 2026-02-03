import logging

logging.basicConfig(filename='workflow.log', level=logging.INFO, format='%(asctime)s %(message)s')

rule extract_umis:
    input:
        r1="{pipeline}/tmp/1_primary/a_trimmed/{sample}_r1.fq.gz",
        r2="{pipeline}/tmp/1_primary/a_trimmed/{sample}_r2.fq.gz"
    output:
        bam=temp("{pipeline}/tmp/1_primary/b_umitagged/{sample}.bam")
    threads: 3
    resources:
        mem_mb=2 * 1024,
        time="01:00:00"
    conda:
        "../../../../../envs/main.yaml"
    log:
        "{pipeline}/logs/extract_umis/{sample}.log"
    params:
        tagsize=3,
        trimsize=3
    shell:
        r"""
        echo "[$(date)] Starting extract_umis for {wildcards.sample}" > {log}
        fgbio -Xmx2g --async-io=true --compression=1 FastqToBam \
            --input {input.r1} {input.r2} \
            --output {output.bam} \
            --read-structures {params.tagsize}M{params.trimsize}S+T {params.tagsize}M{params.trimsize}S+T \
            --sample {wildcards.sample} \
            --library {wildcards.sample} \
        >> {log} 2>&1
        echo "[$(date)] Finished extract_umis for {wildcards.sample}" >> {log}
        """
