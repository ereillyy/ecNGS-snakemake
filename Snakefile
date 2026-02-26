PIPELINE = config["pipeline"]
STATES = config["states"]
FASTQ_DIRS = config["fastq_dirs"]

# Build sample mapping from import sources
SAMPLE_SOURCES = {}  # new_name -> (path, original_name)
ALL_SAMPLES = []
ALL_MATCHED_NORMALS = []

for source_name, source_info in config["import_sources"].items():
    source_path = source_info["path"]
    for original_name, new_name in source_info["samples"].items():
        SAMPLE_SOURCES[new_name] = (source_path, original_name)
        ALL_SAMPLES.append(new_name)

SAMPLE_NORMAL_PAIRS = config["sample_normal_pairs"]
SAMPLES = sorted(set([pair["sample"] for pair in SAMPLE_NORMAL_PAIRS]))
MATCHED_NORMALS = sorted(set([pair["normal"] for pair in SAMPLE_NORMAL_PAIRS]))

def get_sample_normal_pairs():
    return [(pair["sample"], pair["normal"]) for pair in SAMPLE_NORMAL_PAIRS]

def get_sample_source(wildcards):
    sample = wildcards.sample if hasattr(wildcards, 'sample') else wildcards.mn
    return SAMPLE_SOURCES[sample]

include: f"pipelines/{PIPELINE}.smk"

print(f"Using pipeline: {PIPELINE}")
print(f"Samples: {SAMPLES}")
print(f"Matched normals: {MATCHED_NORMALS}")