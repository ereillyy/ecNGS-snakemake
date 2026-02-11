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
        snv_vcf="{pipeline}/tmp/1_primary/e_call/{trimF}r{trimR}/{sample}_snv.vcf",
        indel_vcf="{pipeline}/tmp/1_primary/e_call/{trimF}r{trimR}/{sample}_indel.vcf",
    wildcard_constraints:
        trimF=r"\d+",
        trimR=r"\d+"
    threads: 8
    resources:
        mem_mb=70 * 1024,
        time="12:00:00"
    benchmark:
        "{pipeline}/benchmarks/call/{trimF}r{trimR}/{sample}.tsv"
    conda:
        "../../../../../envs/main.yaml"
    log:
        "{pipeline}/logs/call/{trimF}r{trimR}/{sample}.log"
    shell:
        r"""
        echo "[$(date)] Starting call for {wildcards.sample} with trimF={wildcards.trimF} trimR={wildcards.trimR}" > {log}
        python ~/DupCaller/src/DupCaller.py call \
                --bam {input.sample_bam} \
                --normalBam {input.normal_bam} \
                --reference {input.ref_genome} \
                --threads {threads} \
                --noise {input.noise} \
                --germline {input.snp} \
                --trimF {wildcards.trimF} \
                --trimR {wildcards.trimR} \
                --output {wildcards.pipeline}/tmp/1_primary/e_call/{wildcards.trimF}r{wildcards.trimR}/{wildcards.sample} \
        >> {log} 2>&1
        echo "[$(date)] Finished call for {wildcards.sample}" >> {log}
        """


