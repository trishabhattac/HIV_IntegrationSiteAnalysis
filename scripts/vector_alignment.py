from pathlib import Path
import mappy as mp
import matplotlib.pyplot as plt

# Directories
PROJECT_ROOT = Path(__file__).resolve().parent.parent
VECTOR_FASTA = PROJECT_ROOT / "files" / "input_files" / "pcR_Blunt_II_TOPO.fasta"
FILTERED_DIR = PROJECT_ROOT / "files" / "output_files" / "filtered_reads"
NON_ALIGNED_DIR = PROJECT_ROOT / "files" / "output_files" / "non_aligned_inserts"

def extract_inserts_and_collect_map(fastq_path, out_fasta_path, aligner):
    """
    1. Extracts non-aligning (Human + HIV) inserts to FASTA.
    2. Records read lengths and read-centric alignment boundaries (q_st to q_en).
    """
    insert_count = 0
    read_maps = []  # Stores (read_id, read_len, [(q_st, q_en), ...])

    with open(out_fasta_path, "w") as out_f:
        for read_id, seq, qual in mp.fastx_read(str(fastq_path)):
            read_len = len(seq)
            alignments = list(aligner.map(seq))

            if not alignments:
                # 100% Non-aligned read (pure insert)
                out_f.write(f">{read_id}_full_insert\n{seq}\n")
                insert_count += 1
                read_maps.append((read_id, read_len, []))
            else:
                alignments.sort(key=lambda x: x.q_st)
                
                # Store read query coordinates (q_st, q_en) on the read itself
                aligned_blocks = [(hit.q_st, hit.q_en) for hit in alignments]
                read_maps.append((read_id, read_len, aligned_blocks))

                # Extract non-aligned segments
                last_end = 0
                for idx, hit in enumerate(alignments):
                    if hit.q_st > last_end + 50:  # Minimum insert length cutoff
                        insert_seq = seq[last_end:hit.q_st]
                        out_f.write(f">{read_id}_insert_{idx}\n{insert_seq}\n")
                        insert_count += 1
                    last_end = max(last_end, hit.q_en)

                if read_len > last_end + 50:
                    insert_seq = seq[last_end:]
                    out_f.write(f">{read_id}_insert_tail\n{insert_seq}\n")
                    insert_count += 1

    return insert_count, read_maps


def plot_read_centric_map(read_maps, out_plot_path, max_reads_to_show=40):
    """
    Generates a read-centric sequence map:
    - Grey bar: Entire read span (Highlighting non-aligned insert DNA)
    - Blue overlay: Region aligning to vector backbone
    """
    if not read_maps:
        print("  -> No reads to plot.")
        return

    # Sample top reads for legible figure size
    sampled_reads = read_maps[:max_reads_to_show]

    fig, ax = plt.subplots(figsize=(12, 7))

    for y_idx, (read_id, read_len, vector_blocks) in enumerate(sampled_reads):
        # 1. Base grey line showing the full read (non-aligned region)
        ax.plot([0, read_len], [y_idx, y_idx], color="#B0BEC5", linewidth=5, zorder=1)

        # 2. Blue overlays showing aligned vector backbone segments on the read
        for q_st, q_en in vector_blocks:
            ax.plot([q_st, q_en], [y_idx, y_idx], color="#1E88E5", linewidth=5, zorder=2)

    # Custom legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color="#B0BEC5", lw=4, label="Non-aligned Insert (Human/HIV)"),
        Line2D([0], [0], color="#1E88E5", lw=4, label="Vector Backbone Alignment")
    ]
    ax.legend(handles=legend_elements, loc="upper right", frameon=True)

    ax.set_title(f"Read-Centric Alignment Map (First {len(sampled_reads)} Reads)", fontsize=13, fontweight='bold')
    ax.set_xlabel("Read Length Coordinates (bp)", fontsize=11)
    ax.set_ylabel("Individual Reads", fontsize=11)
    ax.grid(axis='x', linestyle='--', alpha=0.5)

    plt.tight_layout()
    plt.savefig(out_plot_path, dpi=300)
    plt.close()


def main():
    if not VECTOR_FASTA.is_file():
        raise FileNotFoundError(f"\n[ERROR] Vector FASTA not found at:\n{VECTOR_FASTA}\n")

    NON_ALIGNED_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Loading vector backbone: {VECTOR_FASTA.name}...")
    aligner = mp.Aligner(str(VECTOR_FASTA), preset="map-ont")
    if not aligner:
        raise ValueError("Failed to load vector backbone file.")

    first_seq_name = aligner.seq_names[0]
    vector_seq = aligner.seq(first_seq_name)
    print(f"Loaded vector backbone: '{first_seq_name}' ({len(vector_seq):,} bp)\n")

    fastq_files = sorted(FILTERED_DIR.glob("*.fastq"))
    if not fastq_files:
        print(f"Warning: No .fastq files found in {FILTERED_DIR}")
        return

    for fastq_file in fastq_files:
        print(f"Processing {fastq_file.name}...")
        base_name = fastq_file.stem

        out_fasta = NON_ALIGNED_DIR / f"{base_name}_non_aligned_inserts.fasta"
        out_plot = NON_ALIGNED_DIR / f"{base_name}_read_alignment_map.png"

        num_inserts, read_maps = extract_inserts_and_collect_map(fastq_file, out_fasta, aligner)
        print(f"  -> Extracted {num_inserts} non-aligned insert segments.")

        plot_read_centric_map(read_maps, out_plot)
        print(f"  -> Saved read alignment map to {out_plot.name}\n")


if __name__ == "__main__":
    main()