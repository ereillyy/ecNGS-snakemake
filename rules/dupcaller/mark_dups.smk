
rule mark_dups:
    input:
        bam="{pipeline}/tmp/1_primary/c_aligned/{sample}.bam", 
        ref_genome=config["ref"]
    output:
        bam="{pipeline}/tmp/1_primary/d_markdup/{sample}_mkdp.bam",
        dedup_metrics="{pipeline}/qc/opt_dup/{sample}_mkdp_metrics.txt"
    threads: 4
    resources:
        mem_mb=150 * 1024,
        time="05:00:00"
    container:
        config["containers"]["gatk"]
    log:
        "{pipeline}/logs/mark_dups/{sample}.log"
    shell:
        r"""
        echo "[$(date)] Starting mark_dups for {wildcards.sample}" > {log}
        gatk MarkDuplicates \
                -I {input.bam} \
                -O {output.bam} \
                -M {output.dedup_metrics} \
                --READ_NAME_REGEX "(?:.*:)?([0-9]+)[^:]*:([0-9]+)[^:]*:([0-9]+)[^:]*$" \
                --DUPLEX_UMI \
                --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
                --TAGGING_POLICY OpticalOnly \
                --BARCODE_TAG DB \
        >> {log} 2>&1
        echo "[$(date)] Finished mark_dups for {wildcards.sample}" >> {log}
        """


rule mark_dups_index:
    input:
        bam="{pipeline}/tmp/1_primary/d_markdup/{sample}_mkdp.bam"
    output:
        bam_bai="{pipeline}/tmp/1_primary/d_markdup/{sample}_mkdp.bam.bai"
    threads: 4
    resources:
        mem_mb=20 * 1024,
        time="02:00:00"
    conda:
        "../../../../../envs/main.yaml"
    log:
        "{pipeline}/logs/mark_dups_index/{sample}.log"
    shell:
        r"""
        echo "[$(date)] Starting mark_dups_index for {wildcards.sample}" > {log}
             samtools index \
                {input.bam} \
                -@ {threads} \
        >> {log} 2>&1
        echo "[$(date)] Finished mark_dups_index for {wildcards.sample}" >> {log}
        """

