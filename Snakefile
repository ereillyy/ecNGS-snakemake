PIPELINE = config["pipeline"]
# T2_PATH = config["t2_path"]
SAMPLES = list(config["samples"].keys())
MATCHED_NORMALS = config["matched_normals"]
FASTQ_DIRS = config["fastq_dirs"]
STATES = config["states"]
T2_PATH = config["tier2_path"]

def get_normal(wildcards):
    return config["samples"][wildcards.sample]["normal"]

include: f"pipelines/{PIPELINE}.smk"

print(f"Using pipeline: {PIPELINE}")
print(f"Samples: {SAMPLES}")
