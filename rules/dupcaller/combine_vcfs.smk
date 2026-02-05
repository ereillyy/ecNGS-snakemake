rule combine_vcfs_dist2500:
    input:
        snv="{pipeline}/tmp/1_primary/e_call/default/dist2500/{sample}_snv.vcf",
        indel="{pipeline}/tmp/1_primary/e_call/default/dist2500/{sample}_indel.vcf"
    output:
        "{pipeline}/vcf/default/dist2500/{sample}.vcf"
    threads: 1
    resources:
        mem_mb=10 * 1024,
        time="01:00:00"
    localrule: True
    conda:
        "../../../../../envs/main.yaml"
    log:
        "{pipeline}/logs/combine_vcfs/dist2500/{sample}.log"
    shell:
        r"""
        echo "[$(date)] Starting combine_vcfs_dist2500 for {wildcards.sample}" > {log}
        python ../../src/combine_dupcaller_vcfs.py \
            --sample {wildcards.sample} \
            --indir {wildcards.pipeline}/tmp/1_primary/e_call/dist2500 \
            --outdir {wildcards.pipeline}/vcf/dist2500 \
        >> {log} 2>&1
        echo "[$(date)] Finished combine_vcfs_dist2500 for {wildcards.sample}" >> {log}
        """

rule combine_vcfs_no_optical:
    input:
        snv="{pipeline}/tmp/1_primary/e_call/default/no_optical/{sample}_snv.vcf",
        indel="{pipeline}/tmp/1_primary/e_call/default/no_optical/{sample}_indel.vcf"
    output:
        "{pipeline}/vcf/default/no_optical/{sample}.vcf"
    threads: 1
    resources:
        mem_mb=10 * 1024,
        time="01:00:00"
    localrule: True
    conda:
        "../../../../../envs/main.yaml"
    log:
        "{pipeline}/logs/combine_vcfs/no_optical/{sample}.log"
    shell:
        r"""
        echo "[$(date)] Starting combine_vcfs_no_optical for {wildcards.sample}" > {log}
        python ../../src/combine_dupcaller_vcfs.py \
            --sample {wildcards.sample} \
            --indir {wildcards.pipeline}/tmp/1_primary/e_call/no_optical \
            --outdir {wildcards.pipeline}/vcf/no_optical \
        >> {log} 2>&1
        echo "[$(date)] Finished combine_vcfs_no_optical for {wildcards.sample}" >> {log}
        """

rule combine_vcfs_dist100:
    input:
        snv="{pipeline}/tmp/1_primary/e_call/default/dist100/{sample}_snv.vcf",
        indel="{pipeline}/tmp/1_primary/e_call/default/dist100/{sample}_indel.vcf"
    output:
        "{pipeline}/vcf/default/dist100/{sample}.vcf"
    threads: 1
    resources:
        mem_mb=10 * 1024,
        time="01:00:00"
    localrule: True
    conda:
        "../../../../../envs/main.yaml"
    log:
        "{pipeline}/logs/combine_vcfs/dist100/{sample}.log"
    shell:
        r"""
        echo "[$(date)] Starting combine_vcfs_dist100 for {wildcards.sample}" > {log}
        python ../../src/combine_dupcaller_vcfs.py \
            --sample {wildcards.sample} \
            --indir {wildcards.pipeline}/tmp/1_primary/e_call/dist100 \
            --outdir {wildcards.pipeline}/vcf/dist100 \
        >> {log} 2>&1
        echo "[$(date)] Finished combine_vcfs_dist100 for {wildcards.sample}" >> {log}
        """
