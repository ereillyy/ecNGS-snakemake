rule filter_noise:
    input:
        vcf="{pipeline}/tmp/3_vcf/a_consqual/{sample}_cons.vcf",
        noisemask=config["noisemask"],
        chr_sizes=config["chr_sizes"]
    output:
        vcf=temp("{pipeline}/tmp/3_vcf/b_noise/{sample}_cons_noise.vcf")
    threads: 4
    resources:
        mem_mb=50 * 1024,
        time="10:00:00"
    conda:
        "../../../../../envs/main.yaml"
    log:
        "{pipeline}/logs/filter_noise/{sample}.log"
    shell:
        r"""
        echo "[$(date)] Starting filter_noise for {wildcards.sample}" > {log}
        bedtools intersect \
            -a {input.vcf} \
            -b {input.noisemask} \
            -v \
            -header \
            -sorted \
            -g {input.chr_sizes} \
        > {output.vcf} 2>> {log}
        echo "[$(date)] Finished filter_noise for {wildcards.sample}" >> {log}
        """
