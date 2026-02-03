rule filter_mn:
    input:
        vcf="{pipeline}/tmp/3_vcf/b_noise/{sample}_cons_noise.vcf",
        mn_vcf=lambda wc: f"{wc.pipeline}/tmp/2_mn/e_filt_vcf/{config['samples'][wc.sample]['normal']}_filtered.vcf.gz",
        chr_sizes=config["chr_sizes"]
    output:
        vcf=temp("{pipeline}/tmp/3_vcf/c_mn/{sample}_cons_noise_mn.vcf")
    threads: 2
    resources:
        mem_mb=50 * 1024,
        time="10:00:00"
    conda:
        "../../../../../envs/main.yaml"
    log:
        "{pipeline}/logs/filter_mn/{sample}.log"
    shell:
        r"""
        echo "[$(date)] Starting filter_mn for {wildcards.sample}" > {log}
        bedtools intersect \
            -a {input.vcf} \
            -b {input.mn_vcf} \
            -wa \
            -header \
            -sorted \
            -g {input.chr_sizes} \
        | bcftools norm \
            --rm-dup both \
            --output {output.vcf} \
        >> {log} 2>&1
        echo "[$(date)] Finished filter_mn for {wildcards.sample}" >> {log}
        """
