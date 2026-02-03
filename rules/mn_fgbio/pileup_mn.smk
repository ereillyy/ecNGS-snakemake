rule pileup_mn:
    input:
        bam="{pipeline}/tmp/2_mn/c_dedup/{sample}.bam",
        bam_bai="{pipeline}/tmp/2_mn/c_dedup/{sample}.bam.bai",
        ref=config['ref']
    output:
        vcf=temp("{pipeline}/tmp/2_mn/d_unfilt_vcf/{sample}_raw.vcf.gz"),
        tbi=temp("{pipeline}/tmp/2_mn/d_unfilt_vcf/{sample}_raw.vcf.gz.tbi")
    params:
        chr_sizes=config["chr_sizes"]
    threads: 4
    resources:
        mem_mb=50 * 1024,
        time="1-00:00:00"
    conda:
        "../../../../../envs/main.yaml"
    log:
        "{pipeline}/logs/pileup_mn/{sample}.log"
    shell:
        r"""
        echo "[pileup_mn] Starting for {wildcards.sample}" >> {log}
        bcftools mpileup \
            -f {input.ref} \
            -R {params.chr_sizes} \
            -a FORMAT/AD,FORMAT/DP,INFO/AD \
            --max-depth 100000 \
            --threads {threads} \
            --output-type u \
            {input.bam} \
        | bcftools call \
            -mv \
            --threads {threads} \
            --output-type z \
            -o {output.vcf} \
            2>> {log}
        bcftools index -t {output.vcf} 2>> {log}
        echo "[pileup_mn] Finished for {wildcards.sample}" >> {log}
        """
