# HIV-1 Integration Site Analysis

Computational analysis of HIV-1 integration-site positioning relative to transcription start sites (TSSs), with supporting analysis of Oxford Nanopore sequencing data from archived integration-site profiling material.

This repository contains the analysis pipeline developed for an MSc research project examining whether **wild-type (WT) and N74D HIV-1 differ in their genomic proximity to transcription start sites**, and whether this relationship changes following **digoxin treatment**.

## Project overview

The project combines two complementary components:

1. **Integration-site analysis**  
   Published HIV-1 integration-site coordinates from Zhyvoloup *et al.* (2017) were processed and mapped to the nearest annotated human TSS. TSS-distance distributions were compared between WT-Dig0, WT-Dig400, N74D-Dig0 and N74D-Dig400.

2. **ONT sequencing analysis**  
   Oxford Nanopore sequencing reads generated from selected cloned products were assessed for read-length distribution and alignment to the expected pCR™-Blunt II-TOPO™ vector backbone.

The main computational analysis includes **448,282 unique integration sites** across the four experimental groups.

## Biological question

HIV-1 integration is influenced by interactions between the viral capsid and host-cell factors involved in nuclear trafficking and chromatin targeting. The N74D capsid mutation has previously been associated with altered integration-site preference.

This project asks:

> **Does the N74D capsid mutation alter the distance of HIV-1 integration sites from transcription start sites, and is this relationship modified by digoxin treatment?**

Rather than considering only whether an integration site lies within a fixed genomic window, the analysis also treats nearest-TSS distance as a continuous genomic feature.

## Data source

Integration-site coordinates were obtained from Supplementary File S5 of:

> Zhyvoloup A, et al. (2017). *Digoxin reveals a functional connection between HIV-1 integration preference and T-cell activation.* PLOS Pathogens 13(7): e1006460.  
> https://doi.org/10.1371/journal.ppat.1006460

The supplementary workbook contains four experimental groups:

- `AllUIS2014-WT-Dig0`
- `WT-Dig400`
- `N74D-Dig0`
- `N74D-Dig400`

Each condition contains three experimental replicates.

## TSS annotation

The published integration-site coordinates are based on the **hg19/GRCh37** genome assembly. To maintain coordinate compatibility, transcript annotations were obtained from **GENCODE v19**.

Transcript TSS positions were defined in a strand-aware manner:

- `+` strand transcripts: annotated transcript start
- `-` strand transcripts: annotated transcript end

One TSS entry was retained per annotated transcript.

## Computational workflow

The TSS analysis is implemented as a sequential Python workflow.

### 1. Integration-site data preparation

The four worksheets from the published Excel file are imported and combined into a single analysis table.

Processing includes:

- retention of viral genotype, treatment and replicate information
- chromosome-name standardisation
- orientation standardisation
- conversion of genomic positions to numeric coordinates
- checks for missing required fields
- identification of duplicate coordinates for quality control
- assignment of an internal `SiteID`

### 2. GENCODE TSS preparation

The GENCODE v19 GTF annotation is processed to create a transcript-level TSS reference containing genomic coordinates and associated gene/transcript metadata.

### 3. Mapping validation

Before running the full dataset, a fixed random subset of 100 integration sites is used to check coordinate compatibility and nearest-TSS assignment.

### 4. Full nearest-TSS mapping

Integration sites are processed chromosome-by-chromosome using **PyRanges**. Each integration coordinate is assigned to its physically nearest annotated transcript TSS without strand restriction.

The absolute distance between each integration site and its nearest TSS is then calculated in base pairs.

### 5. Distance analysis

Nearest-TSS distances are summarised for each viral genotype and treatment group using:

- median
- first quartile (Q1)
- third quartile (Q3)
- interquartile range
- replicate-level medians

Integration sites are also quantified within cumulative distance thresholds of 1, 5, 10, 25, 50 and 100 kb.

For visualisation, these are represented as non-overlapping intervals:

`0–1 kb`, `1–5 kb`, `5–10 kb`, `10–25 kb`, `25–50 kb`, `50–100 kb`, and `>100 kb`.

The **10-kb threshold** is retained as a literature-aligned comparison with the original Zhyvoloup *et al.* analysis, while the continuous-distance analysis is used to describe the broader distribution without relying solely on a fixed cutoff.

## ONT read-processing workflow

ONT FASTQ files are analysed separately from the TSS pipeline.

The sequencing workflow includes:

1. calculation of read-count and read-length summary statistics
2. visualisation of read-length distributions
3. retention of reads meeting a 3.5-kb length threshold
4. alignment of retained reads to the pCR™-Blunt II-TOPO™ vector reference using **mappy/minimap2**
5. extraction of non-vector sequence segments for further inspection
6. generation of read-centric vector-alignment visualisations

The vector backbone was initially inspected using **SnapGene**, and sequence alignments were independently examined in **Geneious** as a confirmation step.

## Repository organisation

The main analysis scripts are organised by workflow:

```text
scripts/
├── tss_processing/
│   ├── integration-site data preparation
│   ├── GENCODE v19 TSS generation
│   ├── mapping validation
│   ├── full nearest-TSS mapping
│   ├── distance summarisation
│   └── figure generation
│
└── ont_processing/
    ├── read-length analysis
    ├── FASTQ length filtering
    └── vector-alignment analysis
```

Processed tables and generated outputs are kept separately from the scripts so that the workflow can be reproduced from the source data.

## Software

The analysis was developed using:

- Python 3.14.6
- pandas 3.0.5
- PyRanges 0.1.4
- mappy / minimap2
- Matplotlib 3.11.1
- SnapGene
- Geneious

## Running the TSS workflow

The scripts in `scripts/tss_processing/` are intended to be run sequentially, beginning with preparation of the published integration-site workbook and the GENCODE annotation, followed by mapping validation, full nearest-TSS assignment, summary analysis and figure generation.

Before running the analysis, ensure that:

- the published S5 `.xlsx` file is available locally
- `gencode.v19.annotation.gtf.gz` is available locally
- input and output paths in the scripts match the local project structure
- the required Python packages are installed

Example package installation:

```bash
pip install pandas pyranges matplotlib mappy openpyxl
```

## Main findings represented by the analysis

The pipeline was developed to quantify differences in TSS proximity between WT and N74D HIV-1.

The analysis showed that:

- WT integration sites were more TSS-proximal than N74D sites under both treatment conditions
- under Dig0, the median nearest-TSS distance was approximately 5.5 kb for WT and 10.6 kb for N74D
- the WT-N74D difference remained present under Dig400 but was smaller
- the same genotype-associated pattern was maintained across the three experimental replicates
- the continuous TSS-distance analysis provides additional resolution beyond classification using a fixed 10-kb threshold

These results are descriptive; the repository does not treat individual integration sites as independent biological replicates for formal genotype-by-treatment inference.

## Interpretation

TSS proximity is used here as a genomic measure of integration-site positioning. It should not be interpreted as direct evidence that the nearest transcript is transcriptionally active, nor as proof of enrichment relative to random genomic integration.

The strongest inference from this analysis is comparative: **WT and N74D HIV-1 show different nearest-TSS distance profiles within the same published experimental framework.**

## Limitations and future analysis

The current workflow uses hg19/GRCh37 with GENCODE v19 to remain compatible with the original integration coordinates. Future work could remap the published sites to **GRCh38** and repeat the analysis using a contemporary GENCODE annotation.

Additional extensions could include:

- matched random integration-site controls
- integration with the original RNA-seq dataset
- gene-body and promoter annotation
- chromatin-accessibility annotation
- speckle-associated and lamina-associated domain analysis
- linkage of integration-site features with virological measurements
- validation in primary CD4+ T-cell datasets

## Reproducibility

The repository is intended to retain the complete computational workflow used for:

- preprocessing published integration-site coordinates
- constructing the hg19 TSS reference
- validating genomic mapping
- assigning nearest TSSs
- calculating TSS-distance summaries
- generating analysis figures
- processing and assessing ONT sequencing reads

Where possible, analysis parameters and thresholds are defined explicitly within the scripts.

## Author

**Trisha Bhattacharyya**  
MSc Applied Biosciences and Biotechnology  
Imperial College London

## Citation

If using the published integration-site dataset, please cite the original study:

> Zhyvoloup A, et al. (2017). Digoxin reveals a functional connection between HIV-1 integration preference and T-cell activation. *PLOS Pathogens*, 13(7), e1006460. https://doi.org/10.1371/journal.ppat.1006460
