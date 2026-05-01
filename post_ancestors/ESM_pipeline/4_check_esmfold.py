# Check results of esmfold runs
# Some proteins failed due to CUDA out of memory error
from pathlib import Path
import pandas as pd
from tqdm import tqdm

# Load data csv with all sequences
data_csv = Path("../data/data.csv")
data_df = pd.read_csv(data_csv, index_col=0)

# Directory with pdb results
pdb_dir = Path("../esm_pdb")


# Check pdbs for each idx
def get_pdb_path(idx, pdb_dir=pdb_dir):
    return pdb_dir / f"{idx}.pdb"


# List of indices without pdbs
failed_ids = []
for idx in tqdm(data_df.index):
    pdb_path = get_pdb_path(idx)
    if not pdb_path.exists():
        failed_ids.append(idx)

print(f"Number of missing pdbs: {len(failed_ids)} / {len(data_df)}")

# Write a fasta file with all the missing idx
# Use index 169 (last batch is 168)
fasta_file = Path("../fasta/170.fasta")
with fasta_file.open("w") as f:
    for idx in failed_ids:
        seq = data_df["SEQ"][idx]
        f.write(f">{idx}\n{seq}\n")
print(f"Missing pdb fastas written to {fasta_file}")
