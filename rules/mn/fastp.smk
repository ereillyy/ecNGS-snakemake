rule mn_fastp:
    input:
        r1=lambda wc: f"{FASTQ_DIRS['unfiltered']}{wc.sample}_r1.fq.gz",
        r2=lambda wc: f"{FASTQ_DIRS['unfiltered']}{wc.sample}_r2.fq.gz"
    output:
        r1=temp("{pipeline}/tmp/2_mn/a_trimmed/{sample}_r1.fq.gz"),
        r2=temp("{pipeline}/tmp/2_mn/a_trimmed/{sample}_r2.fq.gz"),
        html="{pipeline}/qc/mn_fastp/{sample}.html",
        json="{pipeline}/qc/mn_fastp/{sample}.json"
    wildcard_constraints:
        sample="|".join(MATCHED_NORMALS)
    threads: 4
    resources:
        mem_mb=20 * 1024,
        time="01:00:00"
    conda:
        "../../../../../envs/main.yaml"
    log:
        "{pipeline}/logs/mn_fastp/{sample}.log"
    shell:
        r"""
        echo "[$(date)] Starting fastp for matched normal {wildcards.sample}" > {log}
        fastp \
            -i {input.r1} \
            -I {input.r2} \
            -o {output.r1} \
            -O {output.r2} \
            -h {output.html} \
            -j {output.json} \
            --thread {threads} \
            --trim_front1 6 \
            --trim_front2 6 \
            --qualified_quality_phred 20 \
            --unqualified_percent_limit 40 \
            --n_base_limit 5 \
            --length_required 50 \
            --trim_poly_g \
            --detect_adapter_for_pe \
        >> {log} 2>&1
        echo "[$(date)] Finished fastp for matched normal {wildcards.sample}" >> {log}
        """
