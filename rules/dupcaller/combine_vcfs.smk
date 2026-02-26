
rule combine_vcfs:
    input:
        snv_vcf="{pipeline}/tmp/1_primary/e_call/{normal}/{sample}_snv.vcf",
        indel_vcf="{pipeline}/tmp/1_primary/e_call/{normal}/{sample}_indel.vcf",
    output:
        vcf="{pipeline}/vcf/{normal}/{sample}.vcf"
    threads: 16
    resources:
        mem_mb=10,
        time="24:00:00"
    localrule: True
    conda:
        "../../../../../envs/main.yaml"
    log:
        "{pipeline}/logs/combine_vcfs/{normal}/{sample}.log"
    shell:
        r"""
        echo "[$(date)] Starting combine_vcfs for {wildcards.sample} vs {wildcards.normal}" > {log}
        python ../../src/combine_dupcaller_vcfs.py \
            --sample {wildcards.sample} \
            --indir {wildcards.pipeline}/tmp/1_primary/e_call/{wildcards.normal} \
            --outdir {wildcards.pipeline}/vcf/{wildcards.normal} \
        >> {log} 2>&1
        echo "[$(date)] Finished combine_vcfs for {wildcards.sample} vs {wildcards.normal}" >> {log}
        """
