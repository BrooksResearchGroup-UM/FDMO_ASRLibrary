import argparse
import pandas
import os
import pickle
from ete3 import Tree

def no2id(pyasr_output_dir, cluster_id):
    
    # Get the directory one level up from output_dir
    base_dir = os.path.dirname(pyasr_output_dir)  
    
    # Define the paths to the pickle file and the input tree file
    pickle_file = os.path.join(base_dir, f'no2id_{str(cluster_id)}.pickle')
    input_tree_file = os.path.join(pyasr_output_dir, 'tree-to-reconstruct.newick')
    print("Checking file path:", input_tree_file)
    output_tree_file = os.path.join(pyasr_output_dir, 'tree-to-reconstruct-id.newick')
    print("Checking file path:", output_tree_file)

    # Load the number to ID mapping from the pickle file
    with open(pickle_file, 'rb') as f:
        no_to_id_dict = pickle.load(f)

    # Load the tree from the specified file
    t = Tree(input_tree_file)

    # Replace each leaf's name using the loaded dictionary
    for leaf in t.iter_leaves():
        leaf.name = no_to_id_dict[leaf.name]  # Default to original name if no entry is found

    # Save the modified tree to the new file
    t.write(format=2, outfile=output_tree_file)

def extract_tree(filename):
    trigger_phrase = "tree with node labels for Rod Page's TreeView"
    capture = False
    tree_data = ""
    
    with open(filename, 'r') as file:
        for line in file:
            if trigger_phrase in line:
                # Start capturing from the next line
                capture = True
                continue  # Move to the next iteration to skip the trigger line
            if capture:
                if line.strip():  # Check if the line is not just empty space
                    tree_data += line.strip()
                else:
                    # If an empty line is encountered, assume the end of the tree data
                    break

    return tree_data if tree_data else None  # Return None if no tree data was captured

def parse_arguments():
    parser = argparse.ArgumentParser(description='Run post pyasr processing.')
    parser.add_argument('--output_dir', required=True, help='output directory')
    parser.add_argument('--cluster_id', required=True, help='cluster id')
    return parser.parse_args()

def main():
    args = parse_arguments()

    # Convert seq number to id
    no2id(args.output_dir, args.cluster_id)

    # Get the tree with ancestral nodes
    filename = os.path.join(args.output_dir, 'rst')
    tree_data = extract_tree(filename)
    tree_output = os.path.join(args.output_dir, 'tree-to-reconstruct-node.tree')

    if tree_data:
        # Write the extracted tree data to a file
        with open(tree_output, 'w') as out_file:
            out_file.write(tree_data)
    else:
        print("No tree data found in the file.")

    # Load the tree from the specified file
    t_node = Tree(tree_output)

    # Define the paths to the pickle file and the input tree file
    base_dir = os.path.dirname(args.output_dir)  
    pickle_file = os.path.join(base_dir, f'no2id_{str(args.cluster_id)}.pickle')

    # Load the number to ID mapping from the pickle file
    with open(pickle_file, 'rb') as f:
        no_to_id_dict = pickle.load(f)

    # Get the directory one level up from output_dir
    tree_output_id = os.path.join(args.output_dir, 'tree-to-reconstruct-node-id.tree')

    # Replace each leaf's name using the loaded dictionary
    for leaf in t_node.iter_leaves():
        leaf.name = leaf.name.split('_')[1]
        leaf.name = no_to_id_dict[leaf.name]  # Default to original name if no entry is found

    # Save the modified tree to the new file
    t_node.write(format=2, outfile=tree_output_id)

if __name__ == '__main__':
    main()
