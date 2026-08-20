"""
06_plot_tss_results.py

Generate descriptive figures for the HIV-1 integration-site TSS analysis.

Figures:
    1. Median distance to nearest TSS with Q1-Q3 range
    2. Number of integration sites within discrete TSS-distance bins
    3. Median distance to nearest TSS across experimental replicates

Input files:
    files/processed/tss_analysis/tss_distance_summary_by_group.csv
    files/processed/tss_analysis/tss_distance_windows_by_group.csv
    files/processed/tss_analysis/tss_distance_summary_by_replicate.csv

Output:
    files/processed/tss_analysis/figures/
"""

from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


# PROJECT PATHS

PROJECT_ROOT = Path(__file__).resolve().parents[2]

ANALYSIS_DIR = PROJECT_ROOT / "files" / "processed" / "tss_analysis"
FIGURE_DIR = ANALYSIS_DIR / "figures"

FIGURE_DIR.mkdir(parents=True, exist_ok=True)


GROUP_FILE = ANALYSIS_DIR / "tss_distance_summary_by_group.csv"
WINDOW_FILE = ANALYSIS_DIR / "tss_distance_windows_by_group.csv"
REPLICATE_FILE = ANALYSIS_DIR / "tss_distance_summary_by_replicate.csv"

# LOAD DATA
group_df = pd.read_csv(GROUP_FILE)
window_df = pd.read_csv(WINDOW_FILE)
replicate_df = pd.read_csv(REPLICATE_FILE)

print(f"  Group summary:     {len(group_df)} rows")
print(f"  Window summary:    {len(window_df)} rows")
print(f"  Replicate summary: {len(replicate_df)} rows")

group_order = [
    ("WT", "Dig0"),
    ("WT", "Dig400"),
    ("N74D", "Dig0"),
    ("N74D", "Dig400"),
]

group_labels = {
    ("WT", "Dig0"): "WT-Dig0",
    ("WT", "Dig400"): "WT-Dig400",
    ("N74D", "Dig0"): "N74D-Dig0",
    ("N74D", "Dig400"): "N74D-Dig400",
}

# FIGURE 1
# MEDIAN DISTANCE TO NEAREST TSS

figure_1 = []

for virus, treatment in group_order:

    row = group_df[
        (group_df["Virus"] == virus)
        & (group_df["Treatment"] == treatment)
    ]

    if row.empty:
        continue

    row = row.iloc[0]

    figure_1.append(
        {
            "Group": group_labels[(virus, treatment)],
            "Median": row["Median_distance_bp"],
            "Q1": row["Q1_distance_bp"],
            "Q3": row["Q3_distance_bp"],
        }
    )

figure_1_df = pd.DataFrame(figure_1)

x = np.arange(len(figure_1_df))

lower_error = (
    figure_1_df["Median"] - figure_1_df["Q1"]
)

upper_error = (
    figure_1_df["Q3"] - figure_1_df["Median"]
)

fig, ax = plt.subplots(figsize=(9, 6))

ax.errorbar(
    x,
    figure_1_df["Median"],
    yerr=[lower_error, upper_error],
    fmt="o",
    capsize=5,
    markersize=8,
)

ax.set_xticks(x)
ax.set_xticklabels(figure_1_df["Group"])

ax.set_ylabel("Distance to nearest TSS (bp)")
ax.set_xlabel("Experimental group")
ax.set_title("Median distance of HIV-1 integration sites to nearest TSS")

ax.grid(axis="y", alpha=0.25)

plt.tight_layout()

figure_1_path = FIGURE_DIR / "figure_1_median_tss_distance.png"
plt.savefig(figure_1_path, dpi=300, bbox_inches="tight")
plt.close()

# FIGURE 2
# NON-CUMULATIVE DISTANCE BINS


# The original Step 05 output contains cumulative counts:
#
#   Within_1kb
#   Within_5kb
#   Within_10kb
#   Within_25kb
#   Within_50kb
#   Within_100kb
#
# We convert these into non-cumulative distance bins:
#
#   0-1 kb
#   1-5 kb
#   5-10 kb
#   10-25 kb
#   25-50 kb
#   50-100 kb
#   >100 kb


window_bins = [
    ("0–1 kb", "Within_1kb", None),
    ("1–5 kb", "Within_5kb", "Within_1kb"),
    ("5–10 kb", "Within_10kb", "Within_5kb"),
    ("10–25 kb", "Within_25kb", "Within_10kb"),
    ("25–50 kb", "Within_50kb", "Within_25kb"),
    ("50–100 kb", "Within_100kb", "Within_50kb"),
]


figure_2_rows = []

for _, row in window_df.iterrows():

    virus = row["Virus"]
    treatment = row["Treatment"]

    total_sites = row["Number_of_sites"]

    # Calculate each non-cumulative distance bin.
    for label, cumulative_column, previous_column in window_bins:

        if previous_column is None:
            count = row[cumulative_column]
        else:
            count = (
                row[cumulative_column]
                - row[previous_column]
            )

        figure_2_rows.append(
            {
                "Group": group_labels.get(
                    (virus, treatment),
                    f"{virus}-{treatment}",
                ),
                "Virus": virus,
                "Treatment": treatment,
                "Distance_bin": label,
                "Number_of_sites": int(count),
            }
        )

    # Sites beyond 100 kb.
    within_100kb = row["Within_100kb"]

    figure_2_rows.append(
        {
            "Group": group_labels.get(
                (virus, treatment),
                f"{virus}-{treatment}",
            ),
            "Virus": virus,
            "Treatment": treatment,
            "Distance_bin": ">100 kb",
            "Number_of_sites": int(
                total_sites - within_100kb
            ),
        }
    )


figure_2_df = pd.DataFrame(figure_2_rows)

bin_order = [
    "0–1 kb",
    "1–5 kb",
    "5–10 kb",
    "10–25 kb",
    "25–50 kb",
    "50–100 kb",
    ">100 kb",
]

group_order_labels = [
    "WT-Dig0",
    "WT-Dig400",
    "N74D-Dig0",
    "N74D-Dig400",
]


# Arrange the data so that each group has the same bin order.
plot_df = (
    figure_2_df
    .pivot(
        index="Distance_bin",
        columns="Group",
        values="Number_of_sites",
    )
    .reindex(bin_order)
    .reindex(columns=group_order_labels)
)


fig, ax = plt.subplots(figsize=(11, 6.5))

x = np.arange(len(bin_order))

number_of_groups = len(group_order_labels)
bar_width = 0.18

for i, group in enumerate(group_order_labels):

    offset = (
        i - (number_of_groups - 1) / 2
    ) * bar_width

    ax.bar(
        x + offset,
        plot_df[group],
        width=bar_width,
        label=group,
    )


ax.set_xticks(x)
ax.set_xticklabels(bin_order)

ax.set_xlabel("Distance from nearest TSS")
ax.set_ylabel("Number of integration sites")

ax.set_title(
    "Distribution of HIV-1 integration sites by distance from nearest TSS"
)

ax.legend(
    title="Experimental group",
    frameon=False,
)

ax.grid(
    axis="y",
    alpha=0.25,
)

plt.tight_layout()

figure_2_path = (
    FIGURE_DIR
    / "figure_2_tss_distance_distribution.png"
)

plt.savefig(
    figure_2_path,
    dpi=300,
    bbox_inches="tight",
)

plt.close()

# FIGURE 3
# REPLICATE-LEVEL MEDIAN DISTANCES

print("\nGenerating Figure 3...")
print("  Showing replicate-level median TSS distances by experimental group.")


replicate_plot = replicate_df.copy()

replicate_plot["Group"] = (
    replicate_plot["Virus"].astype(str)
    + "-"
    + replicate_plot["Treatment"].astype(str)
)

replicate_plot["Replicate"] = (
    replicate_plot["Replicate"].astype(int)
)


# Experimental groups on the x-axis
x = np.arange(len(group_order_labels))

# Slight horizontal offsets allow the three replicate points
# within each experimental group to be seen separately.
replicate_offsets = {
    1: -0.12,
    2: 0.00,
    3: 0.12,
}


fig, ax = plt.subplots(figsize=(10, 6))


for replicate in [1, 2, 3]:

    x_values = []
    y_values = []

    for group_index, group in enumerate(group_order_labels):

        row = replicate_plot[
            (replicate_plot["Group"] == group)
            & (replicate_plot["Replicate"] == replicate)
        ]

        if row.empty:
            continue

        x_values.append(
            group_index + replicate_offsets[replicate]
        )

        y_values.append(
            row.iloc[0]["Median_distance_bp"]
        )

    ax.scatter(
        x_values,
        y_values,
        s=70,
        label=f"Replicate {replicate}",
        zorder=3,
    )


ax.set_xticks(x)
ax.set_xticklabels(group_order_labels)

ax.set_xlabel("Experimental group")
ax.set_ylabel("Median distance to nearest TSS (bp)")

ax.set_title(
    "Median TSS distance across experimental replicates"
)

ax.legend(
    title="Experimental replicate",
    frameon=False,
)

ax.grid(
    axis="y",
    alpha=0.25,
)

plt.tight_layout()

figure_3_path = (
    FIGURE_DIR
    / "figure_3_replicate_median_tss_distance.png"
)

plt.savefig(
    figure_3_path,
    dpi=300,
    bbox_inches="tight",
)

plt.close()


# SAVE THE NON-CUMULATIVE BIN DATA

# This is useful because Figure 2 is now based on discrete distance
# bins rather than the original cumulative windows.

bin_output_path = (
    ANALYSIS_DIR
    / "tss_distance_bins_by_group.csv"
)

figure_2_df.to_csv(
    bin_output_path,
    index=False,
)

# --------------------------------------------------------
# Figure 4: Distribution of integration sites by
# non-overlapping nearest-TSS distance bins
# --------------------------------------------------------
windows = pd.read_csv(
    PROJECT_ROOT
    / "files"
    / "processed"
    / "tss_analysis"
    / "tss_distance_windows_by_group.csv"
)

windows["0–1 kb"] = windows["Percent_within_1kb"]

windows["1–5 kb"] = (
    windows["Percent_within_5kb"]
    - windows["Percent_within_1kb"]
)

windows["5–10 kb"] = (
    windows["Percent_within_10kb"]
    - windows["Percent_within_5kb"]
)

windows["10–25 kb"] = (
    windows["Percent_within_25kb"]
    - windows["Percent_within_10kb"]
)

windows["25–50 kb"] = (
    windows["Percent_within_50kb"]
    - windows["Percent_within_25kb"]
)

windows["50–100 kb"] = (
    windows["Percent_within_100kb"]
    - windows["Percent_within_50kb"]
)

windows[">100 kb"] = (
    100 - windows["Percent_within_100kb"]
)

distance_bins = [
    "0–1 kb",
    "1–5 kb",
    "5–10 kb",
    "10–25 kb",
    "25–50 kb",
    "50–100 kb",
    ">100 kb",
]

windows["Bin_total"] = windows[distance_bins].sum(axis=1)

print(
    windows[
        ["Virus", "Treatment", "Bin_total"]
    ]
)
print("\nGenerating Figure 4...")
print("  Showing within-group percentage of integration sites by TSS-distance bin.")

import numpy as np


# Define the order in which experimental groups are plotted
group_order = [
    ("WT", "Dig0"),
    ("WT", "Dig400"),
    ("N74D", "Dig0"),
    ("N74D", "Dig400"),
]

x = np.arange(len(distance_bins))
bar_width = 0.20

fig, ax = plt.subplots(figsize=(10, 6))

for i, (virus, treatment) in enumerate(group_order):

    group_row = windows[
        (windows["Virus"] == virus)
        & (windows["Treatment"] == treatment)
    ].iloc[0]

    percentages = [
        group_row[bin_name]
        for bin_name in distance_bins
    ]

    ax.bar(
        x + (i - 1.5) * bar_width,
        percentages,
        width=bar_width,
        label=f"{virus}-{treatment}",
    )

ax.set_xlabel("Distance from nearest TSS")
ax.set_ylabel("Percentage of integration sites (%)")

ax.set_xticks(x)
ax.set_xticklabels(distance_bins)

ax.legend(title="Experimental group")

fig.tight_layout()

figure_4_path = (
    FIGURE_DIR
    / "figure_4_tss_distance_distribution_percent.png"
)

fig.savefig(
    figure_4_path,
    dpi=300,
    bbox_inches="tight",
)

plt.close(fig)

print(f"  Saved: {figure_4_path}")