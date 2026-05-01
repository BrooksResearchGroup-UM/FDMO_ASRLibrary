# Compute plddts of esmfold predictions
from pathlib import Path
import pandas as pd
from tqdm import trange

# Data csv path
data_csv = "../data/data.csv"
data_df = pd.read_csv(data_csv, index_col=0)

# Retrieve average plddt and ptm from slurm logs
MIN_JOB_IDX = 0
MAX_JOB_IDX = 170
slurm_log_dir = Path("../slurm/esmfold")
assert slurm_log_dir.exists()


def read_slurm_log(
    idx,
    slurm_log_dir=slurm_log_dir,
):
    log_file = slurm_log_dir / f"esmfold-{idx}.out"
    assert log_file.exists()

    score_dict = {}
    for line in log_file.open("r").read().splitlines():
        if "Predicted structure" in line:
            line_split = line.split()
            seq_idx = int(line_split[10])
            plddt = float(line_split[15][:-1])
            ptm = float(line_split[17])

            score_dict[seq_idx] = {
                "plddt": plddt,
                "ptm": ptm,
            }
    return score_dict


score_dict = {}
for idx in trange(MIN_JOB_IDX, MAX_JOB_IDX + 1):
    score_dict.update(read_slurm_log(idx))

# Check all ids in score_dict
for seq_idx in data_df.index:
    if seq_idx not in score_dict.keys():
        print(f"WARNING: {seq_idx} has no score!")

# Write dictionary to csv
score_df = pd.DataFrame(score_dict).T
score_df = score_df.sort_index()
print(score_df)

score_csv = "../plddt/scores.csv"
score_df.to_csv(score_csv)
print(f"Wrote plddts and ptms to {score_csv}")
