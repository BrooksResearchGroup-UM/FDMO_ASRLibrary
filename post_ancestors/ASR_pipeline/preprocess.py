import subprocess
import argparse
import pickle
from Bio import SeqIO

def id2seq(input_fasta, intact_input_fasta, intact_output_fasta, output_phylip, id2no_file, no2id_file):
    unique_sequences = {}
    id2no = {}
    no2id = {}
    seq_num = 1

    # Process the first input FASTA file
    for record in SeqIO.parse(input_fasta, "fasta"):
        sequence = str(record.seq)
        if sequence not in unique_sequences:
            new_id = f"seq{seq_num}"
            id2no[str(record.id)] = new_id
            no2id[new_id] = str(record.id)
            record.id = new_id
            unique_sequences[sequence] = record
            seq_num += 1
        else:
            print('Replicate: ', str(record.id))

    # Process the second input FASTA file
    records = list(SeqIO.parse(intact_input_fasta, "fasta"))
    
    filtered_records = []
    for record in records:
        if record.id in id2no:
            record.id = id2no[record.id]
            record.description = ""
            filtered_records.append(record)

    # Write the sequences to a FASTA file
    with open(intact_output_fasta, "w") as output_handle:
        SeqIO.write(filtered_records, output_handle, "fasta")

    # Write the sequences to a PHYLIP file
    with open(output_phylip, "w") as output_handle:
        SeqIO.write(unique_sequences.values(), output_handle, "phylip")

    # Save the id2seq and seq2id mappings to pickle files
    with open(id2no_file, "wb") as handle:
        pickle.dump(id2no, handle)
    
    with open(no2id_file, "wb") as handle:
        pickle.dump(no2id, handle)

def parse_arguments():
    """Parse command line arguments for input and output base directories."""
    parser = argparse.ArgumentParser(description='Process FASTA files and remove duplicates.')
    parser.add_argument('--input_fasta', required=True, help='Path to the first input FASTA file')
    parser.add_argument('--intact_input_fasta', required=True, help='Path to the second input FASTA file')
    parser.add_argument('--intact_output_fasta', required=True, help='Path for the output FASTA file')
    parser.add_argument('--output_phylip', required=True, help='Path for the output PHYLIP file')
    parser.add_argument('--id2no_file', required=True, help='Path for the ID to sequence number mapping file')
    parser.add_argument('--no2id_file', required=True, help='Path for the sequence number to ID mapping file')
    return parser.parse_args()

def main():
    args = parse_arguments()

    # Run the processing function
    id2seq(args.input_fasta, args.intact_input_fasta, args.intact_output_fasta, args.output_phylip, args.id2no_file, args.no2id_file)

if __name__ == '__main__':
    main()
