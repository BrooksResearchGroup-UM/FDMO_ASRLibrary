# Output blocks of fasta files for running esmfold
import pandas as pd
from pathlib import Path

# Path to data csv
data_csv = Path("../data/data.csv")
assert data_csv.exists()
data_df = pd.read_csv(data_csv, index_col=0)

# Directory for output fasta batches
batch_dir = Path("../fasta")
batch_dir.mkdir(exist_ok=True)

# Output fasta files in blocks of 500
BATCH_SIZE = 500

seqs = data_df["SEQ"]
batches = [seqs[i : i + BATCH_SIZE] for i in range(0, len(seqs), BATCH_SIZE)]

print(f"Batch size: {BATCH_SIZE}, Num Batches: {len(batches)}")

# Write batches to batch directory
seq_idx = 0
for i, batch in enumerate(batches):
    fasta_path = batch_dir / f"{i}.fasta"
    with fasta_path.open("w") as f:
        for seq in batch:
            f.write(f">{seq_idx}\n{seq}\n")
            seq_idx += 1
