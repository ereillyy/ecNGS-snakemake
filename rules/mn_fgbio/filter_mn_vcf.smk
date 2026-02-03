rule filter_mn_vcf:
    input:
        vcf="{pipeline}/tmp/2_mn/d_unfilt_vcf/{sample}_raw.vcf.gz",
        tbi="{pipeline}/tmp/2_mn/d_unfilt_vcf/{sample}_raw.vcf.gz.tbi"
    output:
        vcf="{pipeline}/tmp/2_mn/e_filt_vcf/{sample}_filtered.vcf.gz",
        tbi="{pipeline}/tmp/2_mn/e_filt_vcf/{sample}_filtered.vcf.gz.tbi"
    params:
        cov=config["COV"],
        vaf=config["VAF"]
    threads: 1
    resources:
        mem_mb=10 * 1024,
        time="02:00:00"
    conda:
        "../../../../../envs/main.yaml"
    log:
        "{pipeline}/logs/filter_mn_vcf/{sample}.log"
    shell:
        r"""
        echo "[$(date)] Starting filter_mn_vcf for {wildcards.sample}" > {log}
        bcftools filter \
            -i "FORMAT/DP > {params.cov} && (FORMAT/AD[0:1] / FORMAT/DP) < {params.vaf}" \
            --output-type z \
            -o {output.vcf} \
            {input.vcf} \
        >> {log} 2>&1
        bcftools index -t {output.vcf} >> {log} 2>&1
        echo "[$(date)] Finished filter_mn_vcf for {wildcards.sample}" >> {log}
        """
