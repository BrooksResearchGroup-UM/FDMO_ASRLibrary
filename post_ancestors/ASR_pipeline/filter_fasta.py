import argparse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from ete3 import Tree

def filter_nodes(tree_file, min_taxa):
    
    tree = Tree(tree_file)
    
    nodes = []

    # Traverse tree and evaluate each node
    for node in tree.traverse():
        # print(node.name)
        if not node.is_leaf():  # Check if it is an ancestral node
            num_leaves = len(node)  # This uses the len() on node, which returns the number of leaves under it
            if num_leaves >= min_taxa:
                nodes.append(int(node.support))

    return nodes

def filter_fasta(fasta_file, filter_ids):
    output_fasta_file = fasta_file.strip('.fasta') + '_filtered.fasta'
    with open(fasta_file, 'r') as fasta, open(output_fasta_file, 'w') as output:
        for record in SeqIO.parse(fasta, 'fasta'):
            id = record.id
            des = record.description
            node = id.split('_')[2]
            if node[-1] == 'a':
                node = node[:-1]
            else:
                node = node
            prob = float(des.split(': ')[1])
            # print(prob)
            if int(node) in filter_ids and prob >= 0.85:
                SeqIO.write(record, output, 'fasta')

def parse_arguments():
    parser = argparse.ArgumentParser(description='Filter ancestral sequences.')
    parser.add_argument('--input_tree', required=True, help='Path to the input TREE file')
    parser.add_argument('--input_fasta', required=True, help='Path to the input FASTA file')
    return parser.parse_args()

def main():
    args = parse_arguments()

    ancestral_nodes = filter_nodes(args.input_tree, 5)
    # print(ancestral_nodes)
    # print(len(ancestral_nodes))
    filter_fasta(args.input_fasta, ancestral_nodes[1:])

if __name__ == '__main__':
    main()
