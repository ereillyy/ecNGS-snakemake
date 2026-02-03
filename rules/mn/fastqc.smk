rule mn_fastqc:
    input:
        r1=lambda wc: f"{FASTQ_DIRS[wc.state]}{wc.sample}_r1.fq.gz",
        r2=lambda wc: f"{FASTQ_DIRS[wc.state]}{wc.sample}_r2.fq.gz"
    output:
        r1="{pipeline}/qc/fastqc/{state}/{sample}_r1_fastqc.html",
        r2="{pipeline}/qc/fastqc/{state}/{sample}_r2_fastqc.html"
    wildcard_constraints:
        sample="|".join(MATCHED_NORMALS),
        state="unfiltered|mn_filtered"
    threads: 2
    resources:
        mem_mb=5 * 1024,
        time="00:30:00"
    conda:
        "../../../../../envs/qc.yaml"
    log:
        "{pipeline}/logs/mn_fastqc/{state}/{sample}.log"
    shell:
        """
        echo "[$(date)] Starting fastqc for matched normal {wildcards.sample} ({wildcards.state})" > {log}
        fastqc --threads {threads} \
                --outdir {wildcards.pipeline}/qc/fastqc/{wildcards.state} \
                {input.r1} {input.r2} \
        >> {log} 2>&1
        echo "[$(date)] Finished fastqc for matched normal {wildcards.sample} ({wildcards.state})" >> {log}
        """
