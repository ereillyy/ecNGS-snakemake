import logging

logging.basicConfig(filename='workflow.log', level=logging.INFO, format='%(asctime)s %(message)s')



rule align_umi_consensus:
    input:
        cons="{pipeline}/tmp/1_primary/e_cons/{sample}.bam",
        ref=config['ref']
    output:
        bam=temp("{pipeline}/tmp/1_primary/f_consaligned/{sample}.bam")
    threads: 4
    resources:
        mem_mb=40 * 1024,
        time="04:00:00"
    conda:
        "../../../../../envs/main.yaml"
    log:
        "{pipeline}/logs/align_umi_consensus/{sample}.log"
    shell:
        r"""
        echo "[$(date)] Starting align_umi_consensus for {wildcards.sample}" > {log}
        samtools fastq {input.cons} \
        | bwa-mem2 mem \
            -t {threads} \
            -p \
            -C \
            -K 100000000 \
            -Y {input.ref} - \
        | fgbio -Xmx40g --compression 1 --async-io ZipperBams \
            --unmapped {input.cons} \
            --ref {input.ref} \
            --tags-to-reverse Consensus \
            --tags-to-revcomp Consensus \
            --output {output.bam} \
        >> {log} 2>&1
        echo "[$(date)] Finished align_umi_consensus for {wildcards.sample}" >> {log}
        """
