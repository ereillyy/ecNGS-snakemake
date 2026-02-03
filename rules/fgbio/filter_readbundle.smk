import logging

logging.basicConfig(filename='workflow.log', level=logging.INFO, format='%(asctime)s %(message)s')

rule filter_readbundle:
    input:
        bam="{pipeline}/tmp/1_primary/f_consaligned/{sample}.bam",
        ref=config['ref']
    output:
        filtered=temp("{pipeline}/tmp/2_bam_filtered/a_fgbio/{sample}.bam")
    threads: 3
    resources:
        mem_mb=20 * 1024,
        time="02:00:00"
    conda:
        "../../../../../envs/main.yaml"
    log:
        "{pipeline}/logs/filter_readbundle/{sample}.log"
    shell:
        r"""
        echo "[$(date)] Starting filter_readbundle for {wildcards.sample}" > {log}
        fgbio -Xmx20g FilterConsensusReads \
            --input {input.bam} \
            --ref {input.ref} \
            --min-read 4 2 2 \
            --max-read-error-rate 1 \
            --min-base-quality 10 \
            --max-base-error-rate 1 \
            --max-no-call-fraction 1 \
            --require-single-strand-agreement false \
            --output /dev/stdout \
        | samtools sort \
            --threads {threads} \
            -o {output.filtered} \
        >> {log} 2>&1
        echo "[$(date)] Finished filter_readbundle for {wildcards.sample}" >> {log}
        """
