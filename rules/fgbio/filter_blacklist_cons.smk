rule filter_blacklist_cons:
    input:
        bam="{pipeline}/tmp/2_bam_filtered/a_fgbio/{sample}.bam",
        blacklist=config['blacklist']
    output:
        filtered=temp("{pipeline}/tmp/2_bam_filtered/b_blacklist/{sample}.bam")
    threads: 3
    resources:
        mem_mb=5 * 1024,
        time="00:30:00"
    conda:
        "../../../../../envs/main.yaml"
    log:
        "{pipeline}/logs/filter_blacklist_cons/{sample}.log"
    shell:
        r"""
        echo "[filter_blacklist_cons] Starting for {wildcards.sample}" >> {log}
        bedtools intersect \
            -a {input.bam} \
            -b {input.blacklist} \
            -v \
            > {output.filtered} \
            2>> {log}
        echo "[filter_blacklist_cons] Finished for {wildcards.sample}" >> {log}
        """
