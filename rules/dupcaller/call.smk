
rule call_dist2500:
    input:
        sample_bam="{pipeline}/tmp/1_primary/d_markdup/dist2500/{sample}.bam",
        sample_bam_bai="{pipeline}/tmp/1_primary/d_markdup/dist2500/{sample}.bam.bai",
        normal_bam=lambda wc: f"{wc.pipeline}/tmp/2_mn/c_dedup/{get_normal(wc)}.bam",
        normal_bam_bai=lambda wc: f"{wc.pipeline}/tmp/2_mn/c_dedup/{get_normal(wc)}.bam.bai",
        ref_genome=config["ref"],
        noise=config["noisemask"],
        snp=config["snpmask"]
    output:
        snv_vcf="{pipeline}/tmp/1_primary/e_call/default/dist2500/{sample}_snv.vcf",
        indel_vcf="{pipeline}/tmp/1_primary/e_call/default/dist2500/{sample}_indel.vcf"
    wildcard_constraints:
        sample="(?!" + "|".join(MATCHED_NORMALS) + ").*"
    threads: 16
    resources:
        mem_mb=100 * 1024,
        time="24:00:00"
    conda:
        "../../../../../envs/main.yaml"
    log:
        "{pipeline}/logs/call/dist2500/{sample}.log"
    shell:
        r"""
        echo "[$(date)] Starting call_dist2500 for {wildcards.sample}" > {log}
        python ~/DupCaller/src/DupCaller.py call \
                --bam {input.sample_bam} \
                --normalBam {input.normal_bam} \
                --reference {input.ref_genome} \
                --threads {threads} \
                --noise {input.noise} \
                --germline {input.snp} \
                --output {wildcards.pipeline}/tmp/1_primary/e_call/dist2500/{wildcards.sample} \
        >> {log} 2>&1
        echo "[$(date)] Finished call_dist2500 for {wildcards.sample}" >> {log}
        """

rule call_no_optical:
    input:
        sample_bam="{pipeline}/tmp/1_primary/c_aligned/{sample}.bam",
        sample_bam_bai="{pipeline}/tmp/1_primary/c_aligned/{sample}.bam.bai",
        normal_bam=lambda wc: f"{wc.pipeline}/tmp/2_mn/c_dedup/{get_normal(wc)}.bam",
        normal_bam_bai=lambda wc: f"{wc.pipeline}/tmp/2_mn/c_dedup/{get_normal(wc)}.bam.bai",
        ref_genome=config["ref"],
        noise=config["noisemask"],
        snp=config["snpmask"]
    output:
        snv_vcf="{pipeline}/tmp/1_primary/e_call/default/no_optical/{sample}_snv.vcf",
        indel_vcf="{pipeline}/tmp/1_primary/e_call/default/no_optical/{sample}_indel.vcf"
    wildcard_constraints:
        sample="(?!" + "|".join(MATCHED_NORMALS) + ").*"
    threads: 16
    resources:
        mem_mb=100 * 1024,
        time="24:00:00"
    conda:
        "../../../../../envs/main.yaml"
    log:
        "{pipeline}/logs/call/no_optical/{sample}.log"
    shell:
        r"""
        echo "[$(date)] Starting call_no_optical for {wildcards.sample}" > {log}
        python ~/DupCaller/src/DupCaller.py call \
                --bam {input.sample_bam} \
                --normalBam {input.normal_bam} \
                --reference {input.ref_genome} \
                --threads {threads} \
                --noise {input.noise} \
                --germline {input.snp} \
                --output {wildcards.pipeline}/tmp/1_primary/e_call/no_optical/{wildcards.sample} \
        >> {log} 2>&1
        echo "[$(date)] Finished call_no_optical for {wildcards.sample}" >> {log}
        """

rule call_distance_100:
    input:
        sample_bam="{pipeline}/tmp/1_primary/d_markdup/dist100/{sample}.bam",
        sample_bam_bai="{pipeline}/tmp/1_primary/d_markdup/dist100/{sample}.bam.bai",
        normal_bam=lambda wc: f"{wc.pipeline}/tmp/2_mn/c_dedup/{get_normal(wc)}.bam",
        normal_bam_bai=lambda wc: f"{wc.pipeline}/tmp/2_mn/c_dedup/{get_normal(wc)}.bam.bai",
        ref_genome=config["ref"],
        noise=config["noisemask"],
        snp=config["snpmask"]
    output:
        snv_vcf="{pipeline}/tmp/1_primary/e_call/default/dist100/{sample}_snv.vcf",
        indel_vcf="{pipeline}/tmp/1_primary/e_call/default/dist100/{sample}_indel.vcf"
    wildcard_constraints:
        sample="(?!" + "|".join(MATCHED_NORMALS) + ").*"
    threads: 16
    resources:
        mem_mb=100 * 1024,
        time="24:00:00"
    conda:
        "../../../../../envs/main.yaml"
    log:
        "{pipeline}/logs/call/dist100/{sample}.log"
    shell:
        r"""
        echo "[$(date)] Starting call_dist100 for {wildcards.sample}" > {log}
        python ~/DupCaller/src/DupCaller.py call \
                --bam {input.sample_bam} \
                --normalBam {input.normal_bam} \
                --reference {input.ref_genome} \
                --threads {threads} \
                --noise {input.noise} \
                --germline {input.snp} \
                --output {wildcards.pipeline}/tmp/1_primary/e_call/dist100/{wildcards.sample} \
        >> {log} 2>&1
        echo "[$(date)] Finished call_dist100 for {wildcards.sample}" >> {log}
        """
