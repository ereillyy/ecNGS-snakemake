rule mn_align:
    input:
        r1="{pipeline}/tmp/2_mn/a_trimmed/{sample}_r1.fq.gz",
        r2="{pipeline}/tmp/2_mn/a_trimmed/{sample}_r2.fq.gz",
        ref_genome=config["ref"]
    output:
        bam=temp("{pipeline}/tmp/2_mn/b_aligned/{sample}.bam"),
        bam_bai=temp("{pipeline}/tmp/2_mn/b_aligned/{sample}.bam.bai")
    wildcard_constraints:
        sample="|".join(MATCHED_NORMALS)
    threads: 10
    resources:
        mem_mb=50 * 1024,
        time="10:00:00"
    conda:
        "../../../../../envs/main.yaml"
    log:
        "{pipeline}/logs/mn_align/{sample}.log"
    shell:
        r"""
        echo "[$(date)] Starting align for matched normal {wildcards.sample}" > {log}
        bwa-mem2 mem \
                -t {threads} \
                -R "@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}\tPL:ILLUMINA" \
                {input.ref_genome} \
                {input.r1} {input.r2} \
             | samtools sort \
                - \
                -@ {threads} \
                -o {output.bam} && \
             samtools index \
                {output.bam} \
                -@ {threads} \
        >> {log} 2>&1
        echo "[$(date)] Finished align for matched normal {wildcards.sample}" >> {log}
        """
