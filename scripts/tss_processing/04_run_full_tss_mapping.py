"""
Map HIV-1 integration sites to their nearest transcription start site (TSS).

Purpose
------
Annotate every HIV-1 integration site from the 2017 PLOS Pathogens study
with the nearest GENCODE v19 / hg19 transcript TSS.

For each integration site, the script records:
    - original integration-site metadata
    - nearest TSS coordinate
    - transcript and gene identifiers
    - gene name
    - transcript strand
    - genomic distance to the nearest TSS

This script performs annotation and quality control only.
Statistical analysis and visualisation are performed separately.

Input
-----
files/processed/integration_sites_master.csv
files/processed/hg19_gencode_v19_TSS.csv

Output
------
files/processed/integration_sites_with_TSS.csv
"""

from pathlib import Path
import warnings

import pandas as pd
import pyranges as pr

PROJECT_ROOT = Path(__file__).resolve().parents[2]
PROCESSED_DIR = PROJECT_ROOT / "files" / "processed"

INTEGRATION_FILE = PROCESSED_DIR / "integration_sites_master.csv"
TSS_FILE = PROCESSED_DIR / "hg19_gencode_v19_TSS.csv"
OUTPUT_FILE = PROCESSED_DIR / "integration_sites_with_TSS.csv"

# PyRanges emits a warning about strand columns during its internal join.
# We deliberately perform the nearest-TSS search without strand restriction,
# so this warning does not indicate an analysis problem.

warnings.filterwarnings(
    "ignore",
    message="Strand data from other will be added as strand data to self.*"
)

# REQUIRED COLUMNS

INTEGRATION_COLUMNS = [
    "SiteID",
    "Virus",
    "Treatment",
    "Replicate",
    "CloneID",
    "Chr",
    "Position",
    "Orientation",
    "SourceSheet",
]

TSS_COLUMNS = [
    "Chromosome",
    "TSS",
    "Strand",
    "GeneID",
    "TranscriptID",
    "GeneName",
    "TranscriptName",
]


# LOAD INPUT DATA

integration_df = pd.read_csv(INTEGRATION_FILE)
tss_df = pd.read_csv(TSS_FILE)

# VALIDATE INPUTS
missing_integration = [
    column
    for column in INTEGRATION_COLUMNS
    if column not in integration_df.columns
]

missing_tss = [
    column
    for column in TSS_COLUMNS
    if column not in tss_df.columns
]

if missing_integration:
    raise ValueError(
        "Missing integration-site columns: "
        + ", ".join(missing_integration)
    )

if missing_tss:
    raise ValueError(
        "Missing TSS annotation columns: "
        + ", ".join(missing_tss)
    )
else:
    print("\nInput validation: OK")


# PREPARE GENOMIC COORDINATES

"""
PyRanges uses 0-based, half-open genomic intervals.

A single genomic position is therefore represented as:

    Start = position - 1
    End   = position

The original biological coordinate is retained separately in Position/TSS.
"""

integration_df["Position"] = pd.to_numeric(
    integration_df["Position"],
    errors="raise",
)

tss_df["TSS"] = pd.to_numeric(
    tss_df["TSS"],
    errors="raise",
)

integration_df["Start"] = integration_df["Position"] - 1
integration_df["End"] = integration_df["Position"]

tss_df["Start"] = tss_df["TSS"] - 1
tss_df["End"] = tss_df["TSS"]


# CHECK CHROMOSOME COMPATIBILITY
integration_chromosomes = set(integration_df["Chr"])
tss_chromosomes = set(tss_df["Chromosome"])

shared_chromosomes = sorted(
    integration_chromosomes & tss_chromosomes
)

missing_chromosomes = sorted(
    integration_chromosomes - tss_chromosomes
)

print("\nChromosome compatibility:")
print(f"  Integration chromosomes: {len(integration_chromosomes)}")
print(f"  TSS chromosomes:         {len(tss_chromosomes)}")
print(f"  Shared chromosomes:      {len(shared_chromosomes)}")

if missing_chromosomes:
    raise ValueError(
        "Integration-site chromosomes without TSS annotations: "
        + ", ".join(missing_chromosomes)
    )

print("  All integration chromosomes represented in TSS annotation: OK")


# FIND NEAREST TSS
mapping_results = []

for chromosome in shared_chromosomes:

    integration_chr = integration_df[
        integration_df["Chr"] == chromosome
    ].copy()

    tss_chr = tss_df[
        tss_df["Chromosome"] == chromosome
    ].copy()

    # Create PyRanges object for integration sites.
    integration_ranges = pr.PyRanges(
        integration_chr[
            [
                "Chr",
                "Start",
                "End",
                "SiteID",
                "Virus",
                "Treatment",
                "Replicate",
                "CloneID",
                "Position",
                "Orientation",
                "SourceSheet",
            ]
        ].rename(
            columns={"Chr": "Chromosome"}
        )
    )

    # Create PyRanges object for transcript TSSs.
    tss_ranges = pr.PyRanges(
        tss_chr[
            [
                "Chromosome",
                "Start",
                "End",
                "TSS",
                "Strand",
                "GeneID",
                "TranscriptID",
                "GeneName",
                "TranscriptName",
            ]
        ]
    )

    # Find the physically nearest TSS.
    #
    # strandedness=False means that the nearest TSS is selected regardless
    # of whether the transcript is on the + or - strand.
    nearest = integration_ranges.nearest(
        tss_ranges,
        strandedness=False,
        apply_strand_suffix=False,
    )

    result = nearest.df

    # Calculate the genomic distance explicitly.
    result["Distance_to_TSS_bp"] = (
        result["Position"] - result["TSS"]
    ).abs()

    mapping_results.append(result)

    print(
        f"  {chromosome}: {len(result):,} sites mapped"
    )

# COMBINE RESULTS
mapped_df = pd.concat(
    mapping_results,
    ignore_index=True,
)

# select final dataset
final_columns = [
    "SiteID",
    "Virus",
    "Treatment",
    "Replicate",
    "CloneID",
    "Chromosome",
    "Position",
    "Orientation",
    "SourceSheet",
    "TSS",
    "Strand",
    "GeneID",
    "TranscriptID",
    "GeneName",
    "TranscriptName",
    "Distance_to_TSS_bp",
]

missing_output_columns = [
    column
    for column in final_columns
    if column not in mapped_df.columns
]

if missing_output_columns:
    raise ValueError(
        "Expected output columns missing: "
        + ", ".join(missing_output_columns)
    )

mapped_df = mapped_df[final_columns]


# FINAL QUALITY CONTROL

total_sites = len(integration_df)
mapped_sites = len(mapped_df)

missing_tss = mapped_df["TSS"].isna().sum()
missing_distance = mapped_df["Distance_to_TSS_bp"].isna().sum()

if mapped_sites == total_sites:
    print("Mapping completeness:        100.00%")
else:
    percentage = mapped_sites / total_sites * 100
    print(f"Mapping completeness:        {percentage:.2f}%")

print("\nMissing annotations:")
print(f"  TSS:              {missing_tss:,}")
print(f"  Distance to TSS:  {missing_distance:,}")

# save the final dataset
mapped_df.to_csv(
    OUTPUT_FILE,
    index=False,
)

print(
    "\nThe output contains one row per integration site with "
    "its nearest GENCODE v19 transcript TSS."
)