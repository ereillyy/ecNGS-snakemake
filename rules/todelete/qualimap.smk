# Define BAM processing steps that need sorting before qualimap
STEPS_NEED_SORT = ["1c_aligned", "1d_grouped", "1f_consaligned"]

# Mapping step names to BAM paths
STEP_TO_BAM = {
    "1c_aligned": "tmp/1_primary/c_aligned/{sample}.bam",
    "1d_grouped": "tmp/1_primary/d_grouped/{sample}.bam",
    "1f_consaligned": "tmp/1_primary/f_consaligned/{sample}.bam",
    "2a_bam_filtered": "tmp/2_bam_filtered/a_fgbio/{sample}.bam",
    "2b_blacklist": "tmp/2_bam_filtered/b_blacklist/{sample}.bam",
    "2c_params": "tmp/2_bam_filtered/c_params/{sample}.bam"
}

rule qualimap_sample:
    input:
        bam=lambda wc: f"{wc.pipeline}/{STEP_TO_BAM[wc.step].format(sample=wc.sample)}"
    output:
        html="{pipeline}/qc/qualimap/{step}/{sample}/qualimapReport.html",
        dir=directory("{pipeline}/qc/qualimap/{step}/{sample}/")
    wildcard_constraints:
        step="|".join(STEP_TO_BAM.keys())
    threads: 1
    resources:
        mem_mb=lambda wc: 50 * 1024 if wc.step in STEPS_NEED_SORT else 20 * 1024,
        time="12:00:00"
    conda:
        "../../../../../envs/qc.yaml"
    log:
        "{pipeline}/logs/qualimap/{step}/{sample}.log"
    params:
        needs_sort=lambda wc: wc.step in STEPS_NEED_SORT
    shell:
        r"""
        echo "[$(date)] Starting qualimap_{wildcards.step} for {wildcards.sample}" > {log}
        
        if [ "{params.needs_sort}" = "True" ]; then
            sorted_bam=$(mktemp -u).bam
            samtools sort -o $sorted_bam {input.bam} 2>> {log}
            input_bam=$sorted_bam
        else
            input_bam={input.bam}
        fi
        
        java_mem=$((({resources.mem_mb}/1024)/2))
        qualimap bamqc \
            -bam $input_bam \
            -outformat HTML \
            -outdir {output.dir} \
            --java-mem-size=${{java_mem}}G \
        >> {log} 2>&1
        
        if [ "{params.needs_sort}" = "True" ]; then
            rm $sorted_bam
        fi
        
        echo "[$(date)] Finished qualimap_{wildcards.step} for {wildcards.sample}" >> {log}
        """

# Matched normal qualimap rule
rule qualimap_mn:
    input:
        bam="{pipeline}/tmp/2_mn/c_dedup/{mn}.bam"
    output:
        html="{pipeline}/qc/qualimap/mn/c_dedup/{mn}/qualimapReport.html",
        dir=directory("{pipeline}/qc/qualimap/mn/c_dedup/{mn}/")
    wildcard_constraints:
        mn="|".join(MATCHED_NORMALS)
    threads: 1
    resources:
        mem_mb=20 * 1024,
        time="12:00:00"
    conda:
        "../../../../../envs/qc.yaml"
    log:
        "{pipeline}/logs/qualimap/mn_c_dedup/{mn}.log"
    shell:
        r"""
        echo "[$(date)] Starting qualimap_mn_c_dedup for {wildcards.mn}" > {log}
        java_mem=$((({resources.mem_mb}/1024)/2))
        qualimap bamqc \
            -bam {input.bam} \
            -outformat HTML \
            -outdir {output.dir} \
            --java-mem-size=${{java_mem}}G \
        >> {log} 2>&1
        echo "[$(date)] Finished qualimap_mn_c_dedup for {wildcards.mn}" >> {log}
        """
