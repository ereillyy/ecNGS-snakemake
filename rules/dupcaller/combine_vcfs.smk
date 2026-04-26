

rule combine_vcfs:
    input:
        snv_vcf="{pipeline}/tmp/1_primary/e_call/threshold_{snv}_{indel}/{sample}_snv.vcf",
        indel_vcf="{pipeline}/tmp/1_primary/e_call/threshold_{snv}_{indel}/{sample}_indel.vcf",
    output:
        vcf="{pipeline}/vcf/threshold_{snv}_{indel}/{sample}.vcf"
    threads: 16
    resources:
        mem_mb=10,
        runtime=24 * 60
    localrule: True
    conda:
        "../../../../../envs/main.yaml"
    log:
        "{pipeline}/logs/combine_vcfs/threshold_{snv}_{indel}/{sample}.log"
    shell:
        r"""
        echo "[$(date)] Starting combine_vcfs for {wildcards.sample}" > {log}
        python ../../src/combine_dupcaller_vcfs.py \
            --sample {wildcards.sample} \
            --indir {wildcards.pipeline}/tmp/1_primary/e_call/threshold_{snv}_{indel} \
            --outdir {wildcards.pipeline}/vcf/threshold_{snv}_{indel} \
        >> {log} 2>&1
        echo "[$(date)] Finished combine_vcfs for {wildcards.sample}" >> {log}
        """
