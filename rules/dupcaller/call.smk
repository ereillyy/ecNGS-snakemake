
rule call:
    input:
        sample_bam="{pipeline}/tmp/1_primary/d_markdup/{sample}.bam",
        sample_bam_bai="{pipeline}/tmp/1_primary/d_markdup/{sample}.bam.bai",
        normal_bam=lambda wc: f"{wc.pipeline}/tmp/2_mn/c_dedup/{get_normal(wc)}.bam",
        normal_bam_bai=lambda wc: f"{wc.pipeline}/tmp/2_mn/c_dedup/{get_normal(wc)}.bam.bai",
        ref_genome=config["ref"],
        noise=config["noisemask"],
        snp=config["snpmask"]
    output:
        snv_vcf="{pipeline}/tmp/1_primary/e_call/threshold_{snv}_{indel}/{sample}_snv.vcf",
        indel_vcf="{pipeline}/tmp/1_primary/e_call/threshold_{snv}_{indel}/{sample}_indel.vcf",
    threads: 8
    resources:
        mem_mb=70 * 1024,
        runtime=24 * 60
    conda:
        "../../../../../envs/main.yaml"
    log:
        "{pipeline}/logs/call/threshold_{snv}_{indel}/{sample}.log"
    shell:
        r"""
        echo "[$(date)] Starting call for {wildcards.sample}" > {log}
        python ~/DupCaller/src/DupCaller.py call \
                --bam {input.sample_bam} \
                --normalBam {input.normal_bam} \
                --reference {input.ref_genome} \
                --thresholdSnv {wildcards.snv} \
                --thresholdIndel {wildcards.indel} \
                --threads {threads} \
                --noise {input.noise} \
                --germline {input.snp} \
                --output {wildcards.pipeline}/tmp/1_primary/e_call/threshold_{wildcards.snv}_{wildcards.indel}/{wildcards.sample} \
        >> {log} 2>&1
        echo "[$(date)] Finished call for {wildcards.sample}" >> {log}
        """
