        snv_vcf="{pipeline}/tmp/1_primary/e_call/default/{sample}_snv.vcf",


rule combine_vcfs:
    input:
        snv_vcf="{pipeline}/tmp/1_primary/e_call/default/{sample}_snv.vcf",
        indel_vcf="{pipeline}/tmp/1_primary/e_call/default/{sample}_indel.vcf",
    output:
        vcf="{pipeline}/vcf/default/{sample}.vcf"
    threads: 16
    resources:
        mem_mb=10,
        time="24:00:00"
    localrule: True
    conda:
        "../../../../../envs/main.yaml"
    log:
        "{pipeline}/logs/combine_vcfs/{sample}.log"
    shell:
        r"""
        echo "[$(date)] Starting combine_vcfs for {wildcards.sample}" > {log}
        python ../../src/combine_dupcaller_vcfs.py \
            --sample {wildcards.sample} \
            --indir {wildcards.pipeline}/tmp/1_primary/e_call/default \
            --outdir {wildcards.pipeline}/vcf/default \
        >> {log} 2>&1
        echo "[$(date)] Finished combine_vcfs for {wildcards.sample}" >> {log}
        """
