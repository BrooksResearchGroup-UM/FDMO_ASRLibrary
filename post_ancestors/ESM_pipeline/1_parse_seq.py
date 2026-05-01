# Parse ancestral sequences

# Imports

from pathlib import Path

# Path to ancestral sequence reconstruction results for all clusters:

vae_dir = Path("/home/cdchiang/for_azamh/vae")
pfam_dir = vae_dir / "20240324_2_PF01494"
all_cluster_dir = pfam_dir / "cluster_fasta_40"
print(f"Location of clusters: {all_cluster_dir}\n")

# Location of extant and ancestral sequences for each cluster
anc_format = (
    all_cluster_dir / "cluster_{i}/cluster_{i}_pyasr_model2/AncFDMO_ancestors_{i}.fasta"
)
extant_format = all_cluster_dir / "cluster_{i}.fasta"

CLUSTER_IDX_START = 1
CLUSTER_IDX_END = 40


# Helper function to get fasta paths:
def get_cluster_path(
    idx: str,
    cluster_format: str,
):
    cluster_path = Path(str(cluster_format).format(i=idx))
    return cluster_path


# Check all fasta paths exist for all clusters
def validate_clusters():
    # Validate that all extant and ancestor cluster paths exist
    for cluster_idx in range(CLUSTER_IDX_START, CLUSTER_IDX_END + 1):
        cluster_anc_fasta = get_cluster_path(cluster_idx, anc_format)
        cluster_extant_fasta = get_cluster_path(cluster_idx, extant_format)

        assert cluster_anc_fasta.exists()
        assert cluster_extant_fasta.exists()

        if cluster_idx == 1:
            print("Example: Cluster 1")
            print("Anc Fasta:", cluster_anc_fasta)
            print("Extant Fasta:", cluster_extant_fasta, "\n")


validate_clusters()


# Parse fasta entries for each fasta file
def parse_fasta(fasta_path):
    import biotite.sequence.io.fasta as fasta

    fasta_file = fasta.FastaFile.read(fasta_path)

    headers = []
    sequences = []

    for header, sequence in fasta_file.items():
        headers.append(header)
        sequences.append(sequence)

    return headers, sequences


# Create dataframe containing all sequences
def parse_clusters():
    from tqdm import trange
    import pandas as pd

    data = {}

    print("Parsing sequences into dataframe:")

    for cluster_idx in trange(CLUSTER_IDX_START, CLUSTER_IDX_END):
        cluster_anc_fasta = get_cluster_path(cluster_idx, anc_format)
        cluster_extant_fasta = get_cluster_path(cluster_idx, extant_format)

        anc_headers, anc_sequences = parse_fasta(cluster_anc_fasta)
        extant_headers, extant_sequences = parse_fasta(cluster_extant_fasta)

        types = ["anc"] * len(anc_headers) + ["extant"] * len(extant_headers)
        headers = anc_headers + extant_headers
        seqs = anc_sequences + extant_sequences

        for header, seq, seq_type in zip(headers, seqs, types):
            data[len(data)] = {
                "HEADER": header,
                "SEQ": seq,
                "TYPE": seq_type,
                "CLUSTER": cluster_idx,
            }

    df = pd.DataFrame(data).T
    return df


df = parse_clusters()

# Write dataframe to csv
data_dir = Path("../data").resolve()
data_dir.mkdir(exist_ok=True)
csv_path = data_dir / "data.csv"
df.to_csv(csv_path)
print(f"Wrote dataframe as csv to {csv_path}")
