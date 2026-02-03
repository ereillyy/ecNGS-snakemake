
rule align:
    input:
        r1_umi="{pipeline}/tmp/1_primary/b_umi/{sample}_1.fastq",
        r2_umi="{pipeline}/tmp/1_primary/b_umi/{sample}_2.fastq",
        ref_genome=config["ref"]
    output:
        bam=temp("{pipeline}/tmp/1_primary/c_aligned/{sample}.bam"), 
        bam_bai=temp("{pipeline}/tmp/1_primary/c_aligned/{sample}.bam.bai")
    threads: 10
    resources:
        mem_mb=50 * 1024,
        time="10:00:00"
    conda:
        "../../../../../envs/main.yaml"
    log:
        "{pipeline}/logs/align/{sample}.log"
    shell:
        r"""
        echo "[$(date)] Starting align for {wildcards.sample}" > {log}
        bwa-mem2 mem \
                -C \
                -t {threads} \
                -R "@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}\tPL:ILLUMINA" \
                {input.ref_genome} \
                {input.r1_umi} {input.r2_umi} \
             | samtools sort \
                - \
                -@ {threads} \
                -o {output.bam} && \
             samtools index \
                {output.bam} \
                -@ {threads} \
        >> {log} 2>&1
        echo "[$(date)] Finished align for {wildcards.sample}" >> {log}
        """
