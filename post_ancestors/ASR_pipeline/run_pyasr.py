import phylopandas as pd
import dendropy as d
import pyasr
import argparse
import pandas
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def extract_gamma_shape_parameter(output_file):
    with open(output_file, 'r') as file:
        lines = file.readlines()

    # Initialize a variable to hold the gamma shape parameter
    gamma_shape_parameter = None

    # Iterate through each line in the file
    for line in lines:
        # Check if the line contains the gamma shape parameter
        if "Gamma shape parameter:" in line:
            # Split the line on whitespace and extract the value
            parts = line.split()
            gamma_shape_parameter = parts[-1]  # The last part is the value
            break  # No need to continue through the rest of the lines

    # Return the gamma shape parameter, if found
    return gamma_shape_parameter

def run_pyasr(input_fasta, input_tree, input_stats, output_dir, model, cluster_id):
    # Use phylopandas to read a set of ancestors.
    df_seqs = pd.read_fasta(input_fasta)

    # Use dendropy to read in tree.
    tree = d.Tree.get(path=input_tree, schema='newick')

    # Get the gamma shape parameter
    alpha = extract_gamma_shape_parameter(input_stats)

    # Reconstruct nodes in tree.
    df_seqs, df_anc, tree = pyasr.reconstruct(df_seqs, tree, working_dir=output_dir, alpha=alpha, model=model)

    # Write out ancestor dataframe to a CSV file.
    df_anc.to_csv(f'{output_dir}/AncFDMO_ancestors_{str(cluster_id)}.csv')

    # Convert the DataFrame to a list of SeqRecord objects
    records = []
    for _, row in df_anc.iterrows():
        id = row['id']
        seq = Seq(row['ml_sequence'].replace('-', ''))
        alt_seq = Seq(row['alt_sequence'].replace('-', ''))
        ml_pro = format(row['ml_posterior'], '.3f')
        alt_pro = format(row['alt_posterior'], '.3f')
        
        # Create SeqRecord objects for the sequence and the alternative sequence
        record = SeqRecord(seq, id='AncFDMO' + '_' + str(cluster_id) + '_' + str(id), description=f'posterior probability: {ml_pro}')
        alt_record = SeqRecord(alt_seq, id='AncFDMO' + '_' + str(cluster_id) + '_' + str(id) + 'a', description=f'posterior probability: {alt_pro}')
        
        # Add the SeqRecords to the list
        records.append(record)
        records.append(alt_record)

    # Write the SeqRecords to a FASTA file
    with open(f'{output_dir}/AncFDMO_ancestors_{str(cluster_id)}.fasta', 'w') as file:
        SeqIO.write(records, file, 'fasta')

def parse_arguments():
    parser = argparse.ArgumentParser(description='Run pyasr.')
    parser.add_argument('--input_fasta', required=True, help='Path to the input FASTA file')
    parser.add_argument('--input_tree', required=True, help='Path to the input TREE file')
    parser.add_argument('--input_stats', required=True, help='Path to the input STATS file')
    parser.add_argument('--output_dir', required=True, help='output directory')
    parser.add_argument('--model', required=True, help='model 2 or 3')
    parser.add_argument('--cluster_id', required=True, help='cluster id')
    return parser.parse_args()

def main():
    args = parse_arguments()

    # Run pyasr
    run_pyasr(args.input_fasta, args.input_tree, args.input_stats, args.output_dir, args.model, args.cluster_id)

if __name__ == '__main__':
    main()
