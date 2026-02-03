rule mn:


rule sample:


rule all:
    input:
        expand("{pipeline}/vcf/{sample}.vcf",
        sample=SAMPLES, pipeline=PIPELINE),
        expand("{pipeline}/qc/fastqc/{state}/{sample}_r1_fastqc.html", 
        sample=SAMPLES, state=STATES, pipeline=PIPELINE),
        expand("{pipeline}/tmp/2_mn/e_filt_vcf/{mn}_filtered.vcf.gz",
        mn=MATCHED_NORMALS, pipeline=PIPELINE),
        expand("{pipeline}/qc/fastqc/{state}/{mn}_r1_fastqc.html",
        mn=MATCHED_NORMALS, state=["unfiltered", "mn_filtered"], pipeline=PIPELINE)


# shared rules
include: "../rules/t2_to_t1.smk"
include: "../rules/fastqc.smk"

# sample rules
include: "../rules/fastp.smk"
include: "../rules/fgbio/extract_umis.smk"
include: "../rules/fgbio/align_umi_tagged.smk"
include: "../rules/fgbio/group_umi.smk"
include: "../rules/fgbio/call_umi_consensus.smk"
include: "../rules/fgbio/align_umi_consensus.smk"

include: "../rules/fgbio/filter_readbundle.smk" # outputs error-corrected bam
include: "../rules/fgbio/filter_blacklist_cons.smk"
include: "../rules/fgbio/filter_read_cons.smk"

include: "../rules/fgbio/get_ec_variants.smk"
include: "../rules/fgbio/filter_noise.smk"
include: "../rules/fgbio/filter_mn.smk"
include: "../rules/fgbio/filter_snp.smk"
include: "../rules/fgbio/filter_mutationsonly.smk"

# matched normal rules
include: "../rules/mn/fastqc.smk"
include: "../rules/mn/fastp.smk"
include: "../rules/mn/align.smk"
include: "../rules/mn/remove_dups.smk"
include: "../rules/mn_fgbio/pileup_mn.smk"
include: "../rules/mn_fgbio/filter_mn_vcf.smk"

