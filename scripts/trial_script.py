import mappy as mp

fastq_path = "files/output_files/filtered_reads/77_3_filtered.fastq"

for read_id, seq, qual in mp.fastx_read(fastq_path):
    print(f"Read ID: {read_id}")
    print(f"First 200 bp:\n{seq[:200]}")
    break