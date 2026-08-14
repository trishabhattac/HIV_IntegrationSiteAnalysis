"""
All processing steps have been carried out using S5 supplementary data from Malament et al. (2017), PLoS Pathogens containing 448,282 HIV-1 integration sites. 
The data were downloaded from the PLoS Pathogens website and saved as an Excel file (ppat.1006460.s012.xlsx) in the 'files/input_files' directory of this project.

 Purpose:
   1. Import the four experimental groups from S5.
   2. Combine them into a single analysis dataframe.
   3. Add experimental metadata.
   4. Perform basic quality-control checks.
   5. Save the QC results for reproducibility.

The TSS analysis itself is performed in a later stage.
"""
from pathlib import Path
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]

INPUT_FILE = (
    PROJECT_ROOT/ "files"/ "input_files"/ "ppat.1006460.s012.xlsx")

PROCESSED_DIR = PROJECT_ROOT / "files" / "processed"
QC_DIR = PROJECT_ROOT / "files" / "qc"

PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
QC_DIR.mkdir(parents=True, exist_ok=True)


# 2. EXPERIMENTAL GROUP DEFINITIONS

""" S5 contains four sheets corresponding to the two HIV-1 variants (WT and N74D) under two treatment conditions:
    1. Dig0   = vehicle/control condition
    2. Dig400 = 400 nM digoxin
    
The original paper describes the Dig0 condition as the DMSO vehicle control for the digoxin treatment. """

SHEETS = {
    "AllUIS2014-WT-Dig0": {
        "Virus": "WT",
        "Treatment": "Dig0",
    },
    "AllUIS2014-WT-Dig400": {
        "Virus": "WT",
        "Treatment": "Dig400",
    },
    "AllUIS2014-N74D-Dig0": {
        "Virus": "N74D",
        "Treatment": "Dig0",
    },
    "AllUIS2014-N74D-Dig400": {
        "Virus": "N74D",
        "Treatment": "Dig400",
    },
}

if not INPUT_FILE.exists():
    raise FileNotFoundError(
        f"\nInput file could not be found:\n{INPUT_FILE}\n\n"
        "Please check the filename and project directory structure."
    )

dataframes = []

for sheet_name, metadata in SHEETS.items():

    df = pd.read_excel(
        INPUT_FILE,
        sheet_name=sheet_name,
    )

    # Add experimental metadata so that all four sheets can be
    # combined into one dataframe while retaining their identity.
    df["Virus"] = metadata["Virus"]
    df["Treatment"] = metadata["Treatment"]
    df["SourceSheet"] = sheet_name

    dataframes.append(df)

    print(
        f"{metadata['Virus']:>4} | "
        f"{metadata['Treatment']:>6} | "
        f"{len(df):>8,} integration sites"
    )

sites = pd.concat(
    dataframes,
    ignore_index=True,
)

sites = sites.rename(
    columns={"Ort": "Orientation"}
)

# Remove accidental whitespace from column names.
sites.columns = sites.columns.str.strip()


# BASIC DATA STANDARDISATION

sites["Chr"] = sites["Chr"].astype(str).str.strip()
sites["Orientation"] = sites["Orientation"].astype(str).str.strip()

sites["Position"] = pd.to_numeric(
    sites["Position"],
    errors="coerce",
)

required_columns = [
    "Replicate",
    "CloneID",
    "Chr",
    "Position",
    "Orientation",
    "Virus",
    "Treatment",
    "SourceSheet",
]

missing_columns = [
    column
    for column in required_columns
    if column not in sites.columns
]

if missing_columns:
    raise ValueError(
        "\nThe following expected columns are missing:\n"
        + "\n".join(f"  - {column}" for column in missing_columns)
    )

# DATASET OVERVIEW
total_sites = len(sites)

print(f"\nTotal integration-site records: {total_sites:,}")

print("\nColumns:")
for column in sites.columns:
    print(f"  - {column}")

# MISSING-VALUE CHECK
missing_values = (
    sites[required_columns]
    .isna()
    .sum()
    .rename("Missing_values")
)

missing_values = missing_values.to_frame()

if missing_values["Missing_values"].sum() == 0:
    print("No missing values detected in the required analysis fields.")
else:
    print(missing_values[missing_values["Missing_values"] > 0])


group_counts = (
    sites
    .groupby(["Virus", "Treatment"], sort=True)
    .size()
    .reset_index(name="Number_of_sites")
)

replicate_counts = (
    sites
    .groupby(
        ["Virus", "Treatment", "Replicate"],
        sort=True,
    )
    .size()
    .reset_index(name="Number_of_sites")
)

chromosome_counts = (
    sites["Chr"]
    .value_counts()
    .sort_index()
    .rename_axis("Chromosome")
    .reset_index(name="Number_of_sites")
)

orientation_counts = (
    sites["Orientation"]
    .value_counts(dropna=False)
    .rename_axis("Orientation")
    .reset_index(name="Number_of_sites")
)

# CLONE ID UNIQUENESS
total_clone_ids = sites["CloneID"].nunique()

duplicate_clone_ids = (
    sites["CloneID"].duplicated(keep=False).sum()
)

print(f"Total records:       {total_sites:,}")
print(f"Unique CloneIDs:     {total_clone_ids:,}")
print(f"Records with a duplicated CloneID: {duplicate_clone_ids:,}")

if total_clone_ids == total_sites:
    print("Result: all CloneIDs are unique.")
else:
    print(
        "Result: duplicated CloneIDs are present and should be "
        "investigated before downstream analysis."
    )


# DUPLICATE GENOMIC COORDINATES

"""A genomic coordinate is defined here by chromosome + position.
We do NOT remove duplicates at this stage.

Multiple records at the same coordinate may represent observations from different experimental groups or replicates and therefore should be retained until the biological meaning of duplicates is established.
"""
duplicate_coordinates = sites.duplicated(
    subset=["Chr", "Position"],
    keep=False,
)

number_duplicate_coordinate_records = (
    duplicate_coordinates.sum()
)

number_unique_coordinates = (
    sites[["Chr", "Position"]]
    .drop_duplicates()
    .shape[0]
)

print(
    f"Unique genomic coordinates: "
    f"{number_unique_coordinates:,}"
)

print(
    f"Records belonging to duplicated coordinates: "
    f"{number_duplicate_coordinate_records:,}"
)

# CREATE A CLEAN MASTER DATASET

# Assign an internal identifier to each record.
# This is simply a tracking identifier and does not replace CloneID.

sites.insert(
    0,
    "SiteID",
    range(1, len(sites) + 1),
)

master_columns = [
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

sites = sites[master_columns]

master_file = (
    PROCESSED_DIR
    / "integration_sites_master.csv"
)

sites.to_csv(
    master_file,
    index=False,
)

# SAVE QUALITY-CONTROL TABLES
group_counts.to_csv(
    QC_DIR / "integration_sites_by_experimental_group.csv", index=False,)

replicate_counts.to_csv(
    QC_DIR / "integration_sites_by_replicate.csv",index=False,)

chromosome_counts.to_csv(
    QC_DIR / "integration_sites_by_chromosome.csv",index=False,)

orientation_counts.to_csv(
    QC_DIR / "integration_sites_by_orientation.csv", index=False,)

missing_values.to_csv(
    QC_DIR / "missing_values.csv")


# CREATE A HUMAN-READABLE QC REPORT
qc_report = QC_DIR / "tss_analysis_initial_QC_report.txt"

with open(qc_report, "w") as report:

    report.write("HIV-1 INTEGRATION SITE ANALYSIS\n")
    report.write("Initial Data Import and Quality Control Report\n")
    report.write("=" * 78 + "\n\n")

    report.write(
        "Source dataset:\n"
        f"  {INPUT_FILE.name}\n\n"
    )

    report.write(
        "Total integration-site records:\n"
        f"  {total_sites:,}\n\n"
    )

    report.write("Required columns:\n")
    for column in required_columns:
        report.write(f"  - {column}\n")

    report.write("\nMissing values:\n")
    report.write(
        missing_values.to_string()
    )

    report.write(
        "\n\nIntegration sites by experimental group:\n"
    )
    report.write(
        group_counts.to_string(index=False)
    )

    report.write(
        "\n\nIntegration sites by replicate:\n"
    )
    report.write(
        replicate_counts.to_string(index=False)
    )

    report.write(
        "\n\nChromosomal distribution:\n"
    )
    report.write(
        chromosome_counts.to_string(index=False)
    )

    report.write(
        "\n\nIntegration-site orientation:\n"
    )
    report.write(
        orientation_counts.to_string(index=False)
    )

    report.write(
        "\n\nCloneID assessment:\n"
        f"  Total records: {total_sites:,}\n"
        f"  Unique CloneIDs: {total_clone_ids:,}\n"
        f"  Records with duplicated CloneID: "
        f"{duplicate_clone_ids:,}\n"
    )

    report.write(
        "\nGenomic-coordinate assessment:\n"
        f"  Unique chromosome-position coordinates: "
        f"{number_unique_coordinates:,}\n"
        f"  Records belonging to duplicated coordinates: "
        f"{number_duplicate_coordinate_records:,}\n"
    )