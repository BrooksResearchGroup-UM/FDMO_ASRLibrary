import subprocess
import argparse
from Bio import SeqIO

# Specify the full paths to the executables
modeltest_path = '/home/cdchiang/modeltest/modeltest-ng-static'
phyml_path = '/home/cdchiang/PhyML3.1/PhyML-3.1_linux64'

def run_modeltest(input_file):
    command = [modeltest_path, '-i', input_file, '-h', 'uigf', '-f', 'ef', '-d', 'aa']
    subprocess.run(command, check=True)

def extract_phyml_command_info(input_file):
    with open(input_file, 'r') as file:
        text = file.read()
    lines = text.split('\n')
    for line in lines:
        if line.startswith('  > phyml'):
            phyml_command = line.strip('> ')
            break
    else:
        return None

    args = phyml_command.split()
    phyml_para = {
        'input_file': args[args.index('-i') + 1],
        'data_type': args[args.index('-d') + 1],
        'model': args[args.index('-m') + 1],
        'frequencies': args[args.index('-f') + 1],
        'proportion_invariant': args[args.index('-v') + 1],
        'gamma_distribution': args[args.index('-a') + 1],
        'number_of_categories': args[args.index('-c') + 1],
        'optimization': args[args.index('-o') + 1],
    }

    return phyml_para

def run_phyml(input_file, model, frequencies, proportion_invariant, gamma_distribution):
    command = [phyml_path, '-i', input_file, '-d', 'aa', '-m', model, '-f', frequencies, '-v', proportion_invariant, '-a', gamma_distribution, '-c', '4', '-o', 'tlr', '-s', 'SPR', '--no_memory_check']
    subprocess.run(command, check=True)

def parse_arguments():
    parser = argparse.ArgumentParser(description='Run ModelTest and PhyML on a FASTA file.')
    parser.add_argument('--input_file', required=True, help='Path to the input FASTA file')
    parser.add_argument('--output_modeltest', required=True, help='Output file for ModelTest results')
    return parser.parse_args()

def main():
    args = parse_arguments()

    # Run ModelTest
    run_modeltest(args.input_file)

    # Extract PhyML command information from ModelTest output
    phyml_params = extract_phyml_command_info(args.output_modeltest)

    if phyml_params:
        # Run PhyML with the extracted parameters
        run_phyml(
            phyml_params['input_file'],
            phyml_params['model'],
            phyml_params['frequencies'],
            phyml_params['proportion_invariant'],
            phyml_params['gamma_distribution']
        )
    else:
        print("Error: PhyML command not found in ModelTest output.")

if __name__ == '__main__':
    main()
