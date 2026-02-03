rule group_umi:
    input:
        bam="{pipeline}/tmp/1_primary/c_aligned/{sample}.bam"
    output:
        grouped=temp("{pipeline}/tmp/1_primary/d_grouped/{sample}.bam"),
        hist="{pipeline}/qc/fgbio_familysize/famsize_{sample}.txt"
    threads: 1
    resources:
        mem_mb=50 * 1024,
        time="04:00:00"
    conda:
        "../../../../../envs/main.yaml"
    log:
        "{pipeline}/logs/group_umi/{sample}.log"
    shell:
        r"""
        echo "[$(date)] Starting group_umi for {wildcards.sample}" > {log}
        fgbio -Xmx50g --compression 1 --async-io GroupReadsByUmi \
            --input {input.bam} \
            --strategy Paired \
            --family-size-histogram {output.hist} \
            --output {output.grouped} \
        >> {log} 2>&1
        echo "[$(date)] Finished group_umi for {wildcards.sample}" >> {log}
        """
