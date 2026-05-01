# Filter sequences based on plddt and ptm cutoffs
from pathlib import Path
import pandas as pd

# Read list of sequences
data_df = pd.read_csv("../data/data.csv", index_col=0)

# Read list of scores
score_df = pd.read_csv("../plddt/scores.csv", index_col=0)

# Merge data frames
data_df = pd.concat([data_df, score_df], axis=1)

# Keep only ancestors
data_df = data_df[data_df["TYPE"] == "anc"]

# Drop nan rows
data_df = data_df.dropna(axis=0, how="any")
print(data_df.columns)
print(data_df.index)
print(data_df.shape)

# Define thresholds to filter sequences
PTM_THRESHOLD = 0.9
PLDDT_THRESHOLD = 90
threshold_df = data_df[
    (data_df["plddt"] >= PLDDT_THRESHOLD) & (data_df["ptm"] >= PTM_THRESHOLD)
]
print(f"{PTM_THRESHOLD=}, {PLDDT_THRESHOLD=}")
print(f"Filtered from {len(data_df)} to {len(threshold_df)} seqs")

# Write filtered sequences to file
threshold_df = threshold_df.sort_index(ascending=True)
screen_dir = Path("../screen")
screen_dir.mkdir(exist_ok=True)
screen_data_file = screen_dir / "screen.csv"
threshold_df.to_csv(screen_data_file)

fasta_file = screen_dir / "screen.fasta"
with fasta_file.open("w") as f:
    for idx in threshold_df.index:
        header = threshold_df["HEADER"][idx]
        seq = threshold_df["SEQ"][idx]
        f.write(">{header}\n{seq}\n")

print(f"Wrote filtered sequences to {screen_data_file} and {fasta_file}")
