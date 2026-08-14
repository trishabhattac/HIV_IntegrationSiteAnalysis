"""GENCODE v19 / hg19 TSS ANNOTATION

Purpose:
Extract transcription start sites (TSSs) from the GENCODE v19 GRCh37/hg19 transcript annotation.
Biological definition:
   For transcripts on the + strand, the TSS is the genomic start.
   For transcripts on the - strand, the TSS is the genomic end.

Output:
 A clean table containing one row per annotated transcript and its TSS.
"""
from pathlib import Path
import pandas as pd 


PROJECT_ROOT = Path(__file__).resolve().parents[1]

GTF_FILE = (
    PROJECT_ROOT
    / "files"/ "input_files"/ "gencode.v19.annotation.gtf.gz")

OUTPUT_DIR = PROJECT_ROOT / "files" / "processed"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


OUTPUT_FILE = OUTPUT_DIR / "hg19_gencode_v19_TSS.csv"

if not GTF_FILE.exists():
    raise FileNotFoundError(
        f"\nGENCODE annotation not found:\n{GTF_FILE}"
    )


GTF_COLUMNS = [
    "Chromosome",
    "Source",
    "Feature",
    "Start",
    "End",
    "Score",
    "Strand",
    "Frame",
    "Attributes",
]

gtf = pd.read_csv(
    GTF_FILE,
    sep="\t",
    comment="#",
    header=None,
    names=GTF_COLUMNS,
    compression="gzip",
)

# KEEP TRANSCRIPT FEATURES
"""
The GTF contains genes, transcripts, exons, CDS features, etc.
For this analysis we only need transcript-level records because each transcript has its own transcription start site. 
"""

transcripts = gtf[
    gtf["Feature"] == "transcript"
].copy()

print(f"Transcript records: {len(transcripts):,}")


# EXTRACT GENE AND TRANSCRIPT IDENTIFIERS

""" 
The GTF stores gene/transcript metadata inside the Attributes column.
We extract the identifiers and gene names for easier interpretation.
"""

transcripts["GeneID"] = transcripts["Attributes"].str.extract(
    r'gene_id "([^"]+)"'
)[0]

transcripts["TranscriptID"] = transcripts["Attributes"].str.extract(
    r'transcript_id "([^"]+)"'
)[0]

transcripts["GeneName"] = transcripts["Attributes"].str.extract(
    r'gene_name "([^"]+)"'
)[0]

transcripts["TranscriptName"] = transcripts["Attributes"].str.extract(
    r'transcript_name "([^"]+)"'
)[0]


# CALCULATE TRANSCRIPTION START SITE

""" Transcription proceeds:

  + strand: left → right
              TSS
               ↓
              start --------------------> end

   - strand: right → left
              end
               ↑
       start <-------------------- TSS

 Therefore, the TSS corresponds to:
   + strand → transcript Start
   - strand → transcript End
"""
transcripts["TSS"] = transcripts["Start"]

minus_strand = transcripts["Strand"] == "-"

transcripts.loc[minus_strand, "TSS"] = (
    transcripts.loc[minus_strand, "End"]
)


# KEEP ONLY RELEVANT INFORMATION

tss = transcripts[
    [
        "Chromosome",
        "TSS",
        "Strand",
        "GeneID",
        "GeneName",
        "TranscriptID",
        "TranscriptName",
    ]
].copy()

print(f"Total transcript TSSs: {len(tss):,}")

print(
    f"Unique chromosomes: "
    f"{tss['Chromosome'].nunique()}"
)

print(
    f"Unique genes represented: "
    f"{tss['GeneID'].nunique():,}"
)

print(
    f"Unique transcripts represented: "
    f"{tss['TranscriptID'].nunique():,}"
)


# Check strand distribution.
strand_counts = (
    tss["Strand"]
    .value_counts()
    .rename_axis("Strand")
    .reset_index(name="Number_of_TSSs")
)

print("\nTSSs by strand:")
print(
    strand_counts.to_string(index=False)
)


# Check for missing values in critical fields.
critical_columns = [
    "Chromosome",
    "TSS",
    "Strand",
    "GeneID",
    "TranscriptID",
]

missing = tss[critical_columns].isna().sum()

print("\nMissing values:")
print(
    missing.to_string()
)


# CHECK CHROMOSOME NAMING

unexpected_chromosomes = sorted(
    set(tss["Chromosome"])
    - {
        f"chr{i}" for i in range(1, 23)
    }
    - {"chrX", "chrY"}
)

if unexpected_chromosomes:
    print(f"Additional chromosomes/scaffolds detected: {unexpected_chromosomes}")
else:
    print(
        "Only chr1-chr22, chrX and chrY are present."
    )


# SAVE CLEAN TSS DATASET

tss.to_csv(
    OUTPUT_FILE,
    index=False,
)

print(f"\nOutput file: {OUTPUT_FILE}")
print("\nThis file contains one TSS record per annotated GENCODE v19 transcript.")
