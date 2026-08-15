# Reads the input files and filters reads based on minimum length threshold. The filtered reads are written to the output directory.
from pathlib import Path

INPUT_DIR = Path("files/input_files")
OUTPUT_DIR = Path("files/output_files/filtered_reads")
MIN_LENGTH = 3500


def filter_fastq(input_path: Path, output_path: Path) -> tuple[int, int]:
    total_reads = 0
    kept_reads = 0

    with input_path.open("r") as infile, output_path.open("w") as outfile:
        while True:
            header = infile.readline()
            if not header:
                break
            sequence = infile.readline()
            plus = infile.readline()
            quality = infile.readline()
            if not quality:
                break

            total_reads += 1
            if len(sequence.strip()) >= MIN_LENGTH:
                outfile.write(header)
                outfile.write(sequence)
                outfile.write(plus)
                outfile.write(quality)
                kept_reads += 1

    return total_reads, kept_reads


def main() -> None:
    # Create the output directory if it does not already exist
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Find all FASTQ files (.fastq and .fq) in the input directory
    fastq_files = sorted(INPUT_DIR.glob("*.fastq")) + sorted(INPUT_DIR.glob("*.fq"))

    filtering_results = []

    for fastq_file in fastq_files:
        output_filename = fastq_file.stem.replace("_raw", "_filtered")
        output_path = OUTPUT_DIR / f"{output_filename}.fastq"

        total_reads, retained_reads = filter_fastq(fastq_file, output_path)

        filtering_results.append(
            (fastq_file.name, total_reads, retained_reads)
        )

    if filtering_results:
        print("Read length filter summary:")
        print("File\tTotal Reads\tRetained Reads\tPercent Retained")

        total_reads_processed = 0
        total_reads_retained = 0

        for file_name, total_reads, retained_reads in filtering_results:
            percent_retained = (
                retained_reads / total_reads * 100
                if total_reads else 0.0
            )
            print(
                f"{file_name}\t"
                f"{total_reads}\t"
                f"{retained_reads}\t"
                f"{percent_retained:.1f}%"
            )

            total_reads_processed += total_reads
            total_reads_retained += retained_reads

        overall_retention = (
            total_reads_retained / total_reads_processed * 100
            if total_reads_processed else 0.0
        )
        print(
            f"TOTAL\t"
            f"{total_reads_processed}\t"
            f"{total_reads_retained}\t"
            f"{overall_retention:.1f}%"
        )

    else:
        print("No FASTQ files found in the input directory.")


if __name__ == "__main__":
    main()



