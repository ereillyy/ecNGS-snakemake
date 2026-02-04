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

# Get samples and matched normals from config
SAMPLES = list(config["samples"].keys())
MATCHED_NORMALS = config["matched_normals"]

def get_normal(wildcards):
    return config["samples"][wildcards.sample]["normal"]

def get_sample_source(wildcards):
    """Returns (path, original_name) for a given sample"""
    sample = wildcards.sample if hasattr(wildcards, 'sample') else wildcards.mn
    return SAMPLE_SOURCES[sample]

include: f"pipelines/{PIPELINE}.smk"

print(f"Using pipeline: {PIPELINE}")
print(f"Samples: {SAMPLES}")
print(f"Matched normals: {MATCHED_NORMALS}")