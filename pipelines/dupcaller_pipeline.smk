# Define all trim combinations to test
TRIM_COMBOS = [
    (0, 0),
    (3, 3),
    (0, 5),
    (5, 0),
    (5, 5),
    (5, 7),
    (7, 5),
    (7, 7),
    (7, 12),
    (12, 7),
    (12, 12),
    (30, 30)
]


rule all:
    input:
        # Test all trim combinations
        expand("{pipeline}/vcf/{trimF}r{trimR}/{sample}.vcf",
            sample=SAMPLES, pipeline=PIPELINE,
            trimF=[tc[0] for tc in TRIM_COMBOS],
            trimR=[tc[1] for tc in TRIM_COMBOS]),
        expand("{pipeline}/qc/fastqc/{state}/{sample}_r1_fastqc.html", 
            sample=SAMPLES, state=STATES, pipeline=PIPELINE),
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
        expand("{pipeline}/vcf/{trimF}r{trimR}/{sample}.vcf",
        sample=SAMPLES, pipeline=PIPELINE,
        trimF=[tc[0] for tc in TRIM_COMBOS],
        trimR=[tc[1] for tc in TRIM_COMBOS]),
#        expand("{pipeline}/qc/fastqc/{state}/{sample}_r1_fastqc.html", 
#        sample=SAMPLES, state=STATES, pipeline=PIPELINE),

# shared rules
include: "../rules/t2_to_t1.smk"
include: "../rules/fastqc.smk"

# sample rules
include: "../rules/fastp.smk"
include: "../rules/dupcaller/align.smk"
include: "../rules/dupcaller/extract_umis.smk"
include: "../rules/dupcaller/mark_dups.smk"
include: "../rules/dupcaller/call.smk"
include: "../rules/dupcaller/combine_vcfs.smk"

# matched normal rules
include: "../rules/mn/fastqc.smk"
include: "../rules/mn/fastp.smk"
include: "../rules/mn/align.smk"
include: "../rules/mn/remove_dups.smk"

