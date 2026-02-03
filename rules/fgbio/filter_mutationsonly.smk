rule filter_mutationsonly:
    input:
        vcf="{pipeline}/tmp/3_vcf/d_snp/{sample}_cons_noise_mn_snp.vcf"
    output:
        vcf="{pipeline}/vcf/{sample}.vcf"
    threads: 1
    resources:
        mem_mb=50 * 1024,
        time="10:00:00"
    conda:
        "../../../../../envs/main.yaml"
    log:
        "{pipeline}/logs/filter_mutationsonly/{sample}.log"
    shell:
        r"""
        echo "[$(date)] Starting filter_mutationsonly for {wildcards.sample}" > {log}
        bcftools view \
            --include 'ALT!="."' \
            --output-type v \
            --output-file {output.vcf} \
            {input.vcf} \
        >> {log} 2>&1
        echo "[$(date)] Finished filter_mutationsonly for {wildcards.sample}" >> {log}
        """
