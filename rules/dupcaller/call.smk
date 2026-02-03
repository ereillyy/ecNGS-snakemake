
rule call:
    input:
        sample_bam="{pipeline}/tmp/1_primary/d_markdup/{sample}_mkdp.bam",
        sample_bam_bai="{pipeline}/tmp/1_primary/d_markdup/{sample}_mkdp.bam.bai",
        normal_bam=lambda wc: f"{wc.pipeline}/tmp/2_mn/c_dedup/{get_normal(wc)}.bam",
        normal_bam_bai=lambda wc: f"{wc.pipeline}/tmp/2_mn/c_dedup/{get_normal(wc)}.bam.bai",
        ref_genome=config["ref"],
        noise=config["noisemask"],
        snp=config["snpmask"]
    output:
        snv_vcf="{pipeline}/tmp/1_primary/e_call/default/{sample}_snv.vcf",
        indel_vcf="{pipeline}/tmp/1_primary/e_call/default/{sample}_indel.vcf",
    threads: 16
    resources:
        mem_mb=100 * 1024,
        time="24:00:00"
    conda:
        "../../../../../envs/main.yaml"
    log:
        "{pipeline}/logs/call/{sample}.log"
    shell:
        r"""
        echo "[$(date)] Starting call for {wildcards.sample}" > {log}
        python ~/DupCaller/src/DupCaller.py call \
                --bam {input.sample_bam} \
                --normalBam {input.normal_bam} \
                --reference {input.ref_genome} \
                --threads {threads} \
                --noise {input.noise} \
                --germline {input.snp} \
                --output {wildcards.pipeline}/tmp/1_primary/e_call/default/{wildcards.sample} \
        >> {log} 2>&1
        echo "[$(date)] Finished call for {wildcards.sample}" >> {log}
        """
