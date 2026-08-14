""""
 Sanity-check the mapping of HIV-1 integration sites to the nearest
transcription start site (TSS) using:
      1. Processed HIV-1 integration-site data
      2. GENCODE v19 hg19 TSS annotation
      
This script uses a small random subset first.
Once the mapping is validated, the same approach will be applied to all integration sites.
"""

from pathlib import Path
import pandas as pd
import pyranges as pr

# reading file paths
PROJECT_ROOT = Path(__file__).resolve().parents[1]

INTEGRATION_FILE = (PROJECT_ROOT/ "files"/ "processed"/ "integration_sites_master.csv")

TSS_FILE = (PROJECT_ROOT/ "files"/ "processed"/ "hg19_gencode_v19_TSS.csv")

OUTPUT_FILE = (PROJECT_ROOT/ "files"/ "processed"/ "tss_mapping_sanity_check.csv")


for file in [INTEGRATION_FILE, TSS_FILE]:

    if not file.exists():
        raise FileNotFoundError(
            f"\nRequired file not found:\n{file}"
        )

integration_sites = pd.read_csv(
    INTEGRATION_FILE
)

# CHECK INTEGRATION COORDINATES
integration_sites["Position"] = pd.to_numeric(
    integration_sites["Position"],errors="coerce",)

invalid_positions = (
    integration_sites["Position"].isna().sum()
)

if invalid_positions > 0:

    raise ValueError(
        f"\n{invalid_positions:,} integration sites have invalid genomic coordinates.")


# LOAD TSS ANNOTATION
tss = pd.read_csv(TSS_FILE)

tss["TSS"] = pd.to_numeric(
    tss["TSS"], errors="coerce",)

invalid_tss = (tss["TSS"].isna().sum())

if invalid_tss > 0:
    raise ValueError(
        f"\n{invalid_tss:,} TSS records have invalid genomic coordinates.")

# CHECK CHROMOSOME COMPATIBILITY
integration_chromosomes = set(integration_sites["Chr"].unique())

tss_chromosomes = set(
    tss["Chromosome"].unique()
)

common_chromosomes = (
    integration_chromosomes
    & tss_chromosomes)

print(
    f"  Integration chromosomes: "
    f"{len(integration_chromosomes)}"
)

print(
    f"  TSS chromosomes: "
    f"{len(tss_chromosomes)}"
)

print(
    f"  Shared chromosomes: "
    f"{len(common_chromosomes)}"
)

print(
    f"  {sorted(common_chromosomes)}"
)


""" SELECT A SMALL RANDOM SUBSET:
We intentionally start with 100 sites.
This lets us verify that the genomic mapping is correct before processing all 448,282 integration sites.
"""

SAMPLE_SIZE = 100

sample = integration_sites.sample(
    n=SAMPLE_SIZE,
    random_state=42,
).copy()

print(
    f"\nRandomly selected {len(sample)} "
    "integration sites for mapping."
)

# CONVERT INTEGRATION SITES TO GENOMIC INTERVALS

"""Each integration site represents a single genomic base.
 For PyRanges we represent this as:

   Start = Position
   End   = Position + 1
"""

sample["Start"] = (
    sample["Position"].astype(int)
)

sample["End"] = (
    sample["Position"].astype(int) + 1
)

integration_ranges = pr.PyRanges(
    sample.rename(
        columns={
            "Chr": "Chromosome"
        }
    )[
        [
            "Chromosome",
            "Start",
            "End",
            "Replicate",
            "CloneID",
            "Orientation",
            "Virus",
            "Treatment",
        ]
    ]
)


#  CONVERT TSSs TO GENOMIC INTERVALS

tss["Start"] = (
    tss["TSS"].astype(int)
)

tss["End"] = (
    tss["TSS"].astype(int) + 1
)

tss_ranges = pr.PyRanges(
    tss[
        [
            "Chromosome",
            "Start",
            "End",
            "Strand",
            "GeneID",
            "GeneName",
            "TranscriptID",
            "TranscriptName",
        ]
    ]
)
# find nearest TSS for each integration site

nearest = integration_ranges.nearest(
    tss_ranges,
    strandedness=False,
)


# Convert the result back into a pandas DataFrame.
result = nearest.df


# CALCULATE / RETAIN DISTANCE

if "Distance" not in result.columns:

    raise ValueError(
        "\nPyRanges did not return a Distance column. "
        "Check the installed PyRanges version."
    )

result["Distance_to_TSS_bp"] = (
    result["Distance"].abs()
)

columns_to_display = [
    "Chromosome",
    "Start",
    "Virus",
    "Treatment",
    "GeneName",
    "TranscriptID",
    "Strand",
    "Distance_to_TSS_bp",
]

available_columns = [
    column
    for column in columns_to_display
    if column in result.columns
]

print(
    result[
        available_columns
    ]
    .head(20)
    .to_string(index=False)
)

print(
    f"Integration sites tested: "
    f"{len(sample):,}"
)

# save results to CSV
result.to_csv(
    OUTPUT_FILE,
    index=False,
)