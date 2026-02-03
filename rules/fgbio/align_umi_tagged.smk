import logging

logging.basicConfig(filename='workflow.log', level=logging.INFO, format='%(asctime)s %(message)s')

rule align_umi_tagged:
    input:
        bam="{pipeline}/tmp/1_primary/b_umitagged/{sample}.bam",
        ref=config['ref']
    output:
        bam=temp("{pipeline}/tmp/1_primary/c_aligned/{sample}.bam")
    threads: 10
    resources:
        mem_mb=50 * 1024,
        time="06:00:00"
    conda:
        "../../../../../envs/main.yaml"
    log:
        "{pipeline}/logs/align_umi_tagged/{sample}.log"
    shell:
        r"""
        echo "[$(date)] Starting align_umi_tagged for {wildcards.sample}" > {log}
        samtools fastq {input.bam} \
        | bwa-mem2 mem \
            -t {threads} \
            -p \
            -C \
            -K 100000000 \
            -Y {input.ref} - \
        | fgbio -Xmx50g --compression 1 --async-io ZipperBams \
            --unmapped {input.bam} \
            --ref {input.ref} \
            --output /dev/stdout \
        | samtools sort \
            - \
            --template-coordinate \
            --threads {threads} \
            -o {output.bam} \
        >> {log} 2>&1
        echo "[$(date)] Finished align_umi_tagged for {wildcards.sample}" >> {log}
        """
