rule get_ec_variants:
    input:
        bam="{pipeline}/tmp/2_bam_filtered/c_params/{sample}.bam",
        ref=config["ref"]
    output:
        vcf=temp("{pipeline}/tmp/3_vcf/a_consqual/{sample}_cons.vcf")
    threads: 3
    resources:
        mem_mb=50 * 1024,
        time="08:00:00"
    conda:
        "../../../../../envs/main.yaml"
    log:
        "{pipeline}/logs/get_ec_variants/{sample}.log"
    shell:
        r"""
        echo "[$(date)] Starting get_ec_variants for {wildcards.sample}" > {log}
        bcftools mpileup {input.bam} \
            --fasta-ref {input.ref} \
            --min-MQ 0 \
            --min-BQ 60 \
            --ignore-RG \
            --max-depth 1000000000 \
            --output-type u \
        | bcftools call \
            -P 0.999 \
            --multiallelic-caller \
            --output-type u \
        | bcftools sort \
            --max-mem 30G \
            --output-type v \
            --temp-dir temp_XXXXXX \
        | bcftools norm \
            --rm-dup both  \
        > {output.vcf} 2>> {log}
        echo "[$(date)] Finished get_ec_variants for {wildcards.sample}" >> {log}
        """
