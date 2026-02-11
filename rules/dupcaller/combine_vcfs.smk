snv_vcf="{pipeline}/tmp/1_primary/e_call/{trimF}r{trimR}/{sample}_snv.vcf",


rule combine_vcfs:
    input:
        snv_vcf="{pipeline}/tmp/1_primary/e_call/{trimF}r{trimR}/{sample}_snv.vcf",
        indel_vcf="{pipeline}/tmp/1_primary/e_call/{trimF}r{trimR}/{sample}_indel.vcf",
    output:
        vcf="{pipeline}/vcf/{trimF}r{trimR}/{sample}.vcf"
    wildcard_constraints:
        trimF=r"\d+",
        trimR=r"\d+"
    threads: 16
    resources:
        mem_mb=10,
        time="24:00:00"
    localrule: True
    conda:
        "../../../../../envs/main.yaml"
    log:
        "{pipeline}/logs/combine_vcfs/{trimF}r{trimR}/{sample}.log"
    shell:
        r"""
        echo "[$(date)] Starting combine_vcfs for {wildcards.sample}" > {log}
        python ../../src/combine_dupcaller_vcfs.py \
            --sample {wildcards.sample} \
            --indir {wildcards.pipeline}/tmp/1_primary/e_call/{wildcards.trimF}r{wildcards.trimR} \
            --outdir {wildcards.pipeline}/vcf/{wildcards.trimF}r{wildcards.trimR} \
        >> {log} 2>&1
        echo "[$(date)] Finished combine_vcfs for {wildcards.sample}" >> {log}
        """
