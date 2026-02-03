rule fastp:
    input:
        r1=lambda wc: f"{FASTQ_DIRS['unfiltered']}{wc.sample}_r1.fq.gz",
        r2=lambda wc: f"{FASTQ_DIRS['unfiltered']}{wc.sample}_r2.fq.gz"
    output:
        r1_trim=temp("{pipeline}/tmp/1_primary/a_trimmed/{sample}_r1.fq.gz"),
        r2_trim=temp("{pipeline}/tmp/1_primary/a_trimmed/{sample}_r2.fq.gz"),
        json="{pipeline}/qc/fastp_reports/{sample}_fastp.json",
        html="{pipeline}/qc/fastp_reports/{sample}_fastp.html"
    threads: 4
    resources:
        mem_mb=20 * 1024,
        time="01:00:00"
    conda:
        "../../../../envs/main.yaml"
    log:
        "{pipeline}/logs/fastp/{sample}.log"
    shell:
        r"""
        echo "[$(date)] Starting fastp for {wildcards.sample}" > {log}
        fastp -i {input.r1} -I {input.r2} \
            -o {output.r1_trim} \
            -O {output.r2_trim} \
            --thread {threads} \
            --disable_quality_filtering \
            --dup_calc_accuracy 4 \
            --trim_poly_g \
            --detect_adapter_for_pe \
            --length_required 15 \
            --json {output.json} \
            --html {output.html} \
        >> {log} 2>&1
        echo "[$(date)] Finished fastp for {wildcards.sample}" >> {log}
        """

