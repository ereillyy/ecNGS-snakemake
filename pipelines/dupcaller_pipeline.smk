rule testimport:
    input:
        expand(f"{FASTQ_DIRS['unfiltered']}{{sample}}_r1.fq.gz",
            sample=SAMPLES)

rule all:
    input:
        # All sample-normal pairs for VCFs
        expand("{pipeline}/vcf/{normal}/{sample}.vcf",
            sample=[pair[0] for pair in get_sample_normal_pairs()],
            normal=[pair[1] for pair in get_sample_normal_pairs()],
            pipeline=PIPELINE),
        #expand("{pipeline}/qc/fastqc/{state}/{sample}_r1_fastqc.html", 
        #    sample=SAMPLES, state=STATES, pipeline=PIPELINE),
        expand("{pipeline}/qc/fastqc/{state}/{mn}_r1_fastqc.html",
            mn=MATCHED_NORMALS, state=["unfiltered", "mn_filtered"], pipeline=PIPELINE),
        expand("{pipeline}/tmp/2_mn/c_dedup/{mn}.bam",
            mn=MATCHED_NORMALS, pipeline=PIPELINE)

rule mn:
    input:
        expand("{pipeline}/tmp/2_mn/c_dedup/{mn}.bam",
        mn=MATCHED_NORMALS, pipeline=PIPELINE),
        expand("{pipeline}/qc/fastqc/{state}/{mn}_r1_fastqc.html",
        mn=MATCHED_NORMALS, state=["unfiltered", "mn_filtered"], pipeline=PIPELINE)

rule sample:
    input:
        expand("{pipeline}/vcf/{normal}/{sample}.vcf",
        sample=[pair[0] for pair in get_sample_normal_pairs()],
        normal=[pair[1] for pair in get_sample_normal_pairs()],
        pipeline=PIPELINE),
        #expand("{pipeline}/qc/fastqc/{state}/{sample}_r1_fastqc.html", 
        #sample=SAMPLES, state=STATES, pipeline=PIPELINE),

# shared rules
#include: "../rules/t2_to_t1.smk"
include: "../rules/fastqc.smk"

# sample rules
#include: "../rules/fastp.smk"
#include: "../rules/dupcaller/extract_umis.smk"
#include: "../rules/dupcaller/align.smk"
#include: "../rules/dupcaller/mark_dups.smk"
include: "../rules/dupcaller/call.smk"
include: "../rules/dupcaller/combine_vcfs.smk"

# matched normal rules
include: "../rules/mn/fastqc.smk"
include: "../rules/mn/fastp.smk"
include: "../rules/mn/align.smk"
include: "../rules/mn/remove_dups.smk"

