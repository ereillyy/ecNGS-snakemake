rule filter_snp:
    input:
        vcf="{pipeline}/tmp/3_vcf/c_mn/{sample}_cons_noise_mn.vcf",
        snpmask=config["snpmask"],
        chr_sizes=config["chr_sizes"]
    output:
        vcf="{pipeline}/tmp/3_vcf/d_snp/{sample}_cons_noise_mn_snp.vcf"
    threads: 4
    resources:
        mem_mb=50 * 1024,
        time="10:00:00"
    conda:
        "../../../../../envs/main.yaml"
    log:
        "{pipeline}/logs/filter_snp/{sample}.log"
    shell:
        r"""
        echo "[$(date)] Starting filter_snp for {wildcards.sample}" > {log}
        bedtools intersect \
            -a {input.vcf} \
            -b {input.snpmask} \
            -v \
            -header \
            -sorted \
            -g {input.chr_sizes} \
        > {output.vcf} 2>> {log}
        echo "[$(date)] Finished filter_snp for {wildcards.sample}" >> {log}
        """
