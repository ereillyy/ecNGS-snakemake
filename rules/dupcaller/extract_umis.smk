rule extract_umis:
    input:
        r1_trim="{pipeline}/tmp/1_primary/a_trimmed/{sample}_r1.fq.gz",
        r2_trim="{pipeline}/tmp/1_primary/a_trimmed/{sample}_r2.fq.gz",
    output:
        r1_umi=temp("{pipeline}/tmp/1_primary/b_umi/{sample}_1.fastq"),
        r2_umi=temp("{pipeline}/tmp/1_primary/b_umi/{sample}_2.fastq")
    threads: 4
    resources:
        mem_mb=10 * 1024,
        time="02:00:00"
    conda:
        "../../../../../envs/main.yaml"
    log:
        "{pipeline}/logs/extract_umis/{sample}.log"
    shell:
        r"""
        echo "[$(date)] Starting extract_umis for {wildcards.sample}" > {log}
        python ~/DupCaller/src/DupCaller.py trim \
                -i {input.r1_trim} \
                -i2 {input.r2_trim} \
                -p NNNXXX \
                -o {wildcards.pipeline}/tmp/1_primary/b_umi/{wildcards.sample} \
        >> {log} 2>&1
        echo "[$(date)] Finished extract_umis for {wildcards.sample}" >> {log}
        """
