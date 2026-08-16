from pathlib import Path

import csv
import matplotlib.pyplot as plt

INPUT_DIR = Path("files/input_files")
OUTPUT_DIR = Path("files/output_files")
FIGURE_DIR = OUTPUT_DIR / "figures"

MIN_LENGTH = 3500

# FASTQ processing

def analyse_fastq(input_path: Path) -> dict:
    """
    Read a FASTQ file and calculate read-length statistics.

    No sequences are modified or filtered during this step.
    """

    read_lengths = []

    with input_path.open("r") as infile:

        while True:
            header = infile.readline()

            if not header:
                break

            sequence = infile.readline()
            plus = infile.readline()
            quality = infile.readline()

            if not quality:
                break

            read_lengths.append(len(sequence.strip()))

    if not read_lengths:
        return {
            "File": input_path.name,
            "Total_reads": 0,
            "Minimum_length_bp": 0,
            "Q1_length_bp": 0,
            "Median_length_bp": 0,
            "Mean_length_bp": 0,
            "Q3_length_bp": 0,
            "Maximum_length_bp": 0,
            "Reads_below_3500bp": 0,
            "Reads_at_or_above_3500bp": 0,
            "Percent_retained_3500bp": 0,
            "Read_lengths": [],
        }

    sorted_lengths = sorted(read_lengths)
    total_reads = len(sorted_lengths)

    # Calculate quartiles using linear interpolation
    def percentile(values, percentile):
        index = (len(values) - 1) * percentile
        lower = int(index)
        upper = min(lower + 1, len(values) - 1)
        fraction = index - lower

        return (
            values[lower]
            + fraction * (values[upper] - values[lower])
        )

    q1 = percentile(sorted_lengths, 0.25)
    median = percentile(sorted_lengths, 0.50)
    q3 = percentile(sorted_lengths, 0.75)

    reads_below_threshold = sum(
        length < MIN_LENGTH for length in read_lengths
    )

    reads_at_or_above_threshold = sum(
        length >= MIN_LENGTH for length in read_lengths
    )

    percent_retained = (
        reads_at_or_above_threshold / total_reads * 100
    )

    return {
        "File": input_path.name,
        "Total_reads": total_reads,
        "Minimum_length_bp": min(sorted_lengths),
        "Q1_length_bp": q1,
        "Median_length_bp": median,
        "Mean_length_bp": sum(sorted_lengths) / total_reads,
        "Q3_length_bp": q3,
        "Maximum_length_bp": max(sorted_lengths),
        "Reads_below_3500bp": reads_below_threshold,
        "Reads_at_or_above_3500bp": reads_at_or_above_threshold,
        "Percent_retained_3500bp": percent_retained,
        "Read_lengths": read_lengths,
    }


def save_summary(results: list[dict], output_path: Path) -> None:

    fieldnames = [
        "File",
        "Total_reads",
        "Minimum_length_bp",
        "Q1_length_bp",
        "Median_length_bp",
        "Mean_length_bp",
        "Q3_length_bp",
        "Maximum_length_bp",
        "Reads_below_3500bp",
        "Reads_at_or_above_3500bp",
        "Percent_retained_3500bp",
    ]

    with output_path.open("w", newline="") as outfile:

        writer = csv.DictWriter(
            outfile,
            fieldnames=fieldnames
        )

        writer.writeheader()

        for result in results:

            row = {
                key: result[key]
                for key in fieldnames
            }

            writer.writerow(row)


# Generate read-length distribution figure

def generate_read_length_figure(results: list[dict], output_path: Path) -> None:

    plt.figure(figsize=(9, 6))

    for result in results:

        lengths = result["Read_lengths"]

        if not lengths:
            continue

        plt.hist(
            lengths,
            bins=40,
            alpha=0.45,
            label=result["File"].replace("_raw.fastq", "")
        )

    # Vector backbone threshold
    plt.axvline(
        MIN_LENGTH,
        linestyle="--",
        linewidth=1.5,
        label="3500 bp threshold"
    )

    plt.xscale("log")

    plt.xlabel("Read length (bp)")
    plt.ylabel("Number of reads")

    plt.title("Read-length distribution of ONT sequencing datasets")

    plt.legend(
        frameon=False
    )

    plt.tight_layout()

    plt.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()


def main() -> None:

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    FIGURE_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    fastq_files = (
        sorted(INPUT_DIR.glob("*.fastq"))
        + sorted(INPUT_DIR.glob("*.fq"))
    )

    if not fastq_files:
        print("No FASTQ files found in the input directory.")
        return

    results = []


    for fastq_file in fastq_files:
        result = analyse_fastq(fastq_file)

        results.append(result)

    # Save CSV
    summary_path = OUTPUT_DIR / "read_length_summary.csv"

    save_summary(
        results,
        summary_path
    )

    # Generate figure
    figure_path = (
        FIGURE_DIR
        / "read_length_distribution.png"
    )

    generate_read_length_figure(
        results,
        figure_path
    )

    for result in results:

        print(
            f"{result['File']}\t"
            f"{result['Total_reads']}\t"
            f"{result['Median_length_bp']:.1f}\t"
            f"{result['Mean_length_bp']:.1f}\t"
            f"{result['Reads_at_or_above_3500bp']}\t"
            f"{result['Percent_retained_3500bp']:.1f}%"
        )


if __name__ == "__main__":
    main()