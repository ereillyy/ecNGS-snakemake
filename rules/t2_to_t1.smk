rule t2_to_t1:
    input:
        r1=lambda wc: f"{get_sample_source(wc)[0]}{get_sample_source(wc)[1]}_r1.fq.gz",
        r2=lambda wc: f"{get_sample_source(wc)[0]}{get_sample_source(wc)[1]}_r2.fq.gz"
    output:
        r1=temp(f"{FASTQ_DIRS['unfiltered']}{{sample}}_r1.fq.gz"),
        r2=temp(f"{FASTQ_DIRS['unfiltered']}{{sample}}_r2.fq.gz")
    threads: 2
    resources:
        io=1
    localrule: True
    log:
        "t2_to_t1_logs/{sample}.log"
    shell:
        r"""
        echo "[$(date)] [t2_to_t1] Starting for {wildcards.sample}" > {log}
        echo "[$(date)] [t2_to_t1] Source: {input.r1}" >> {log}
        rsync -Pah {input.r1} {input.r2} $(dirname {output.r1})/ >> {log} 2>&1
        echo "[$(date)] [t2_to_t1] Finished for {wildcards.sample}" >> {log}
        """