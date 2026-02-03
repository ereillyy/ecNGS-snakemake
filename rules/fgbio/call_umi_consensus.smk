import logging

logging.basicConfig(filename='workflow.log', level=logging.INFO, format='%(asctime)s %(message)s')

rule call_umi_consensus:
    input:
        grouped="{pipeline}/tmp/1_primary/d_grouped/{sample}.bam"
    output:
        cons=temp("{pipeline}/tmp/1_primary/e_cons/{sample}.bam")
    threads: 3
    resources:
        mem_mb=10 * 1024,
        time="02:00:00"
    conda:
        "../../../../../envs/main.yaml"
    log:
        "{pipeline}/logs/call_umi_consensus/{sample}.log"
    shell:
        r"""
        echo "[$(date)] Starting call_umi_consensus for {wildcards.sample}" > {log}
        fgbio -Xmx10g --compression 1 --async-io CallDuplexConsensusReads \
            --input {input.grouped} \
            --consensus-call-overlapping-bases true \
            --error-rate-pre-umi 45 \
            --error-rate-post-umi 40 \
            --min-input-base-quality 10 \
            --trim false \
            --min-reads 1 \
            --output {output.cons} \
            --threads {threads} \
        >> {log} 2>&1
        echo "[$(date)] Finished call_umi_consensus for {wildcards.sample}" >> {log}
        """
