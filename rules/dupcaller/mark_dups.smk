rule mark_dups_dist2500:
    """Mark optical duplicates with 2500 pixel distance (NovaSeq X standard)"""
    input:
        bam="{pipeline}/tmp/1_primary/c_aligned/{sample}.bam"
    output:
        bam="{pipeline}/tmp/1_primary/d_markdup/dist2500/{sample}.bam",
        metrics="{pipeline}/qc/opt_dup/dist2500/{sample}_metrics.txt"
    wildcard_constraints:
        sample="(?!" + "|".join(MATCHED_NORMALS) + ").*"
    threads: 4
    resources:
        mem_mb=150 * 1024,
        time="04:00:00"
    container:
        config["containers"]["gatk"]
    log:
        "{pipeline}/logs/mark_dups/dist2500/{sample}.log"
    shell:
        r"""
        echo "[$(date)] Starting mark_dups_dist2500 for {wildcards.sample}" > {log}
        gatk MarkDuplicates \
            --INPUT {input.bam} \
            --OUTPUT {output.bam} \
            --METRICS_FILE {output.metrics} \
            --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
            --TAGGING_POLICY OpticalOnly \
            --BARCODE_TAG RX \
        >> {log} 2>&1
        echo "[$(date)] Finished mark_dups_dist2500 for {wildcards.sample}" >> {log}
        """

rule index_dist2500:
    """Index BAM after marking duplicates with dist2500"""
    input:
        bam="{pipeline}/tmp/1_primary/d_markdup/dist2500/{sample}.bam"
    output:
        bai="{pipeline}/tmp/1_primary/d_markdup/dist2500/{sample}.bam.bai"
    wildcard_constraints:
        sample="(?!" + "|".join(MATCHED_NORMALS) + ").*"
    threads: 2
    resources:
        mem_mb=4 * 1024,
        time="01:00:00"
    conda:
        "../../../../../envs/main.yaml"
    log:
        "{pipeline}/logs/index/dist2500/{sample}.log"
    shell:
        r"""
        echo "[$(date)] Indexing BAM for {wildcards.sample}" > {log}
        samtools index {input.bam} >> {log} 2>&1
        echo "[$(date)] Finished indexing for {wildcards.sample}" >> {log}
        """

rule mark_dups_dist100:
    """Mark optical duplicates with 100 pixel distance (old unpatterned flowcell default)"""
    input:
        bam="{pipeline}/tmp/1_primary/c_aligned/{sample}.bam"
    output:
        bam="{pipeline}/tmp/1_primary/d_markdup/dist100/{sample}.bam",
        metrics="{pipeline}/qc/opt_dup/dist100/{sample}_metrics.txt"
    wildcard_constraints:
        sample="(?!" + "|".join(MATCHED_NORMALS) + ").*"
    threads: 4
    resources:
        mem_mb=150 * 1024,
        time="04:00:00"
    container:
        config["containers"]["gatk"]
    log:
        "{pipeline}/logs/mark_dups/dist100/{sample}.log"
    shell:
        r"""
        echo "[$(date)] Starting mark_dups_dist100 for {wildcards.sample}" > {log}
        gatk MarkDuplicates \
            --INPUT {input.bam} \
            --OUTPUT {output.bam} \
            --METRICS_FILE {output.metrics} \
            --OPTICAL_DUPLICATE_PIXEL_DISTANCE 100 \
            --TAGGING_POLICY OpticalOnly \
            --BARCODE_TAG RX \
        >> {log} 2>&1
        echo "[$(date)] Finished mark_dups_dist100 for {wildcards.sample}" >> {log}
        """

rule index_dist100:
    """Index BAM after marking duplicates with dist100"""
    input:
        bam="{pipeline}/tmp/1_primary/d_markdup/dist100/{sample}.bam"
    output:
        bai="{pipeline}/tmp/1_primary/d_markdup/dist100/{sample}.bam.bai"
    wildcard_constraints:
        sample="(?!" + "|".join(MATCHED_NORMALS) + ").*"
    threads: 2
    resources:
        mem_mb=4 * 1024,
        time="01:00:00"
    conda:
        "../../../../../envs/main.yaml"
    log:
        "{pipeline}/logs/index/dist100/{sample}.log"
    shell:
        r"""
        echo "[$(date)] Indexing BAM for {wildcards.sample}" > {log}
        samtools index {input.bam} >> {log} 2>&1
        echo "[$(date)] Finished indexing for {wildcards.sample}" >> {log}
        """

