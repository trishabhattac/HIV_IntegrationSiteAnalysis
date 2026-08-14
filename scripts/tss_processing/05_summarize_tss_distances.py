"""
05_summarize_tss_distances.py

Generate descriptive summaries of the distances between HIV-1 integration
sites and their nearest transcription start sites (TSS).

Input:
    files/processed/integration_sites_with_TSS.csv

Outputs:
    files/processed/tss_analysis/
        tss_distance_summary_by_group.csv
        tss_distance_windows_by_group.csv
        tss_distance_summary_by_replicate.csv
        tss_distance_summary.txt

This script performs descriptive analysis only.
No statistical hypothesis testing or figures are generated here.
"""

from pathlib import Path
import pandas as pd


# File paths
PROJECT_ROOT = Path(__file__).resolve().parents[2]

INPUT_FILE = (
    PROJECT_ROOT
    / "files"
    / "processed"
    / "integration_sites_with_TSS.csv"
)

OUTPUT_DIR = (
    PROJECT_ROOT
    / "files"
    / "processed"
    / "tss_analysis"
)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

DISTANCE_WINDOWS = [1_000, 5_000, 10_000, 25_000, 50_000, 100_000]
GROUP_COLUMNS = ["Virus", "Treatment"]

# read the input file
df = pd.read_csv(INPUT_FILE)

required_columns = {
    "Virus",
    "Treatment",
    "Replicate",
    "Chromosome",
    "Distance_to_TSS_bp",
}

missing_columns = required_columns - set(df.columns)

if missing_columns:
    raise ValueError(
        f"Missing required columns: {sorted(missing_columns)}"
    )


# Remove rows without a TSS distance, if any
df = df.dropna(subset=["Distance_to_TSS_bp"]).copy()

df["Distance_to_TSS_bp"] = pd.to_numeric(
    df["Distance_to_TSS_bp"],
    errors="coerce",
)

df = df.dropna(subset=["Distance_to_TSS_bp"]).copy()

total_sites = len(df)
unique_sites = df["SiteID"].nunique() if "SiteID" in df.columns else total_sites
unique_chromosomes = df["Chromosome"].nunique()
unique_genes = df["GeneID"].nunique() if "GeneID" in df.columns else "N/A"
unique_transcripts = (
    df["TranscriptID"].nunique()
    if "TranscriptID" in df.columns
    else "N/A"
)

# 1. Distance summary by experimental group

group_summary = (
    df.groupby(GROUP_COLUMNS)["Distance_to_TSS_bp"]
    .agg(
        Number_of_sites="count",
        Mean_distance_bp="mean",
        Median_distance_bp="median",
        SD_distance_bp="std",
        Q1_distance_bp=lambda x: x.quantile(0.25),
        Q3_distance_bp=lambda x: x.quantile(0.75),
        Minimum_distance_bp="min",
        Maximum_distance_bp="max",
    )
    .reset_index()
)

group_summary["IQR_distance_bp"] = (
    group_summary["Q3_distance_bp"]
    - group_summary["Q1_distance_bp"]
)

# Round values for readability
for column in [
    "Mean_distance_bp",
    "Median_distance_bp",
    "SD_distance_bp",
    "Q1_distance_bp",
    "Q3_distance_bp",
    "IQR_distance_bp",
]:
    group_summary[column] = group_summary[column].round(1)


group_summary.to_csv(
    OUTPUT_DIR / "tss_distance_summary_by_group.csv",
    index=False,
)


# 2. Distance windows by experimental group

window_rows = []

for (virus, treatment), group in df.groupby(GROUP_COLUMNS):

    n_sites = len(group)
    distances = group["Distance_to_TSS_bp"]

    row = {
        "Virus": virus,
        "Treatment": treatment,
        "Number_of_sites": n_sites,
    }

    for window in DISTANCE_WINDOWS:

        count = (distances <= window).sum()
        percentage = (count / n_sites) * 100

        window_label = f"{window // 1000}kb"

        row[f"Within_{window_label}"] = int(count)
        row[f"Percent_within_{window_label}"] = round(
            percentage,
            2,
        )

    window_rows.append(row)


window_summary = pd.DataFrame(window_rows)

window_summary.to_csv(
    OUTPUT_DIR / "tss_distance_windows_by_group.csv",
    index=False,
)


# 3. Replicate-level summary

replicate_summary = (
    df.groupby(["Virus", "Treatment", "Replicate"])
    ["Distance_to_TSS_bp"]
    .agg(
        Number_of_sites="count",
        Median_distance_bp="median",
        Q1_distance_bp=lambda x: x.quantile(0.25),
        Q3_distance_bp=lambda x: x.quantile(0.75),
    )
    .reset_index()
)

for column in [
    "Median_distance_bp",
    "Q1_distance_bp",
    "Q3_distance_bp",
]:
    replicate_summary[column] = replicate_summary[column].round(1)


replicate_summary.to_csv(
    OUTPUT_DIR / "tss_distance_summary_by_replicate.csv",
    index=False,
)


# 4. Create human-readable TXT summary
summary_file = OUTPUT_DIR / "tss_distance_summary.txt"

with open(summary_file, "w", encoding="utf-8") as file:

    file.write("HIV-1 INTEGRATION SITE → TSS DISTANCE ANALYSIS\n")
    file.write("PRELIMINARY SUMMARY\n\n")

    file.write("Input\n")
    file.write("integration_sites_with_TSS.csv\n\n")

    file.write("Total integration sites\n")
    file.write(f"{total_sites:,}\n\n")

    file.write("Dataset\n")
    file.write(f"Unique integration sites: {unique_sites:,}\n")
    file.write(f"Unique chromosomes: {unique_chromosomes:,}\n")
    file.write(f"Unique genes: {unique_genes:,}\n")
    file.write(f"Unique transcripts: {unique_transcripts:,}\n\n")

    # --------------------------------------------------------
    # Group summary
    # --------------------------------------------------------

    file.write("TSS distance by experimental group\n\n")

    file.write(
        "Virus   Treatment    Sites    Mean (bp)    Median (bp)    "
        "SD (bp)    Q1 (bp)    Q3 (bp)    Min (bp)    Max (bp)\n"
    )

    for _, row in group_summary.iterrows():

        file.write(
            f"{row['Virus']:<7}"
            f"{row['Treatment']:<12}"
            f"{int(row['Number_of_sites']):>8}"
            f"{row['Mean_distance_bp']:>13,.1f}"
            f"{row['Median_distance_bp']:>15,.1f}"
            f"{row['SD_distance_bp']:>13,.1f}"
            f"{row['Q1_distance_bp']:>12,.1f}"
            f"{row['Q3_distance_bp']:>12,.1f}"
            f"{row['Minimum_distance_bp']:>12,.0f}"
            f"{row['Maximum_distance_bp']:>12,.0f}\n"
        )

    file.write("\n")

    # --------------------------------------------------------
    # Distance windows
    # --------------------------------------------------------

    file.write("Distance windows\n\n")

    file.write(
        "Virus   Treatment    Sites    "
        "≤1 kb    ≤5 kb    ≤10 kb    ≤25 kb    ≤50 kb    ≤100 kb\n"
    )

    for _, row in window_summary.iterrows():

        file.write(
            f"{row['Virus']:<7}"
            f"{row['Treatment']:<12}"
            f"{int(row['Number_of_sites']):>8}"
            f"{row['Percent_within_1kb']:>9.2f}%"
            f"{row['Percent_within_5kb']:>9.2f}%"
            f"{row['Percent_within_10kb']:>10.2f}%"
            f"{row['Percent_within_25kb']:>10.2f}%"
            f"{row['Percent_within_50kb']:>10.2f}%"
            f"{row['Percent_within_100kb']:>11.2f}%\n"
        )

    file.write("\n")

    # --------------------------------------------------------
    # Replicate summary
    # --------------------------------------------------------

    file.write("Replicate-level medians\n\n")

    file.write(
        "Virus   Treatment    Replicate    Sites    "
        "Median (bp)    Q1 (bp)    Q3 (bp)\n"
    )

    for _, row in replicate_summary.iterrows():

        file.write(
            f"{row['Virus']:<7}"
            f"{row['Treatment']:<12}"
            f"{int(row['Replicate']):>8}"
            f"{int(row['Number_of_sites']):>10}"
            f"{row['Median_distance_bp']:>15,.1f}"
            f"{row['Q1_distance_bp']:>12,.1f}"
            f"{row['Q3_distance_bp']:>12,.1f}\n"
        )


print("DESCRIPTIVE TSS DISTANCE ANALYSIS COMPLETE")