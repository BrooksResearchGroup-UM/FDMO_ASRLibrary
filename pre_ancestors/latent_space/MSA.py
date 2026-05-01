import pickle
import sys
import numpy as np
import os
from sys import exit
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

## change current directory (Change this!!!)
os.chdir("/home/cdchiang/for_azamh/vae")

## define output directory (Change this!!!)
os.mkdir("./data/test")
output_folder = "data/test"

## A MSA file including both the training and testing dataset should be loaded here.
'''The testing dataset including WT and evolved seuqneces was pasted to the end of the training dataset 
and combined into a single file. All sequences were aligned by hmmer-3.1 using a "seed" file of 
protein family 01494 (PF01494). You are able to download the seed file from Pfam database'''

file_name = "./data/PF01494_MSA_4.fasta"

## TropB was selected to be the query sequence for downstream data processing.
query_seq_id = "B8M9J8" #TropB

## read all the sequences into a dictionary
id_list = []
seq_list = []
seq_dict = {}
for record in SeqIO.parse(file_name, "fasta"):
    ID = record.id
    id_list.append(ID)
    seq_list.append(str(record.seq.upper()))
    seq_dict[ID] = str(record.seq.upper())

FMO = []
for idx, (seq, ID) in enumerate(zip(seq_list, id_list)):
    FMO.append(SeqRecord(Seq(seq), id=ID, description=""))
    
## Usually, the description of each seuqnece is super long so I want to remove it from the input and save it to a different fasta file
with open("./{}/PF01494_MSA_5.fasta".format(output_folder), "w") as handle:
    SeqIO.write(FMO, handle, "fasta")
    
## remove gaps in the query sequences
query_seq = seq_dict[query_seq_id] ## with gaps
idx = [ s == "-" or s == "." for s in query_seq]
for k in seq_dict.keys():
    seq_dict[k] = [seq_dict[k][i] for i in range(len(seq_dict[k])) if idx[i] == False]
query_seq = seq_dict[query_seq_id] ## without gaps

## remove sequences with too many gaps (including more than 20% gap)
len_query_seq = len(query_seq)
seq_id = list(seq_dict.keys())
num_gaps = []
removed_seq = []
for k in seq_id:
    num_gaps.append(seq_dict[k].count("-") + seq_dict[k].count("."))
    if seq_dict[k].count("-") + seq_dict[k].count(".") > 0.2*len_query_seq:
        seq_dict.pop(k)
        removed_seq.append(k)
with open("./{}/seq_dict.pkl".format(output_folder), 'wb') as file_handle:
     pickle.dump(seq_dict, file_handle)

## convert aa type into num 0-20
aa = ['R', 'H', 'K',
      'D', 'E',
      'S', 'T', 'N', 'Q',
      'C', 'G', 'P',
      'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W']
aa_index = {}
aa_index['-'] = 0
aa_index['.'] = 0
i = 1
for a in aa:
    aa_index[a] = i
    i += 1
    
## Encode 20 amino acid including gap from 0 to 20
with open("./{}/aa_index.pkl".format(output_folder), 'wb') as file_handle:
    pickle.dump(aa_index, file_handle)
    
seq_msa = []
keys_list = []
for k in seq_dict.keys():
    if seq_dict[k].count('X') > 0 or seq_dict[k].count('Z') or seq_dict[k].count('O')> 0:
        continue    
    seq_msa.append([aa_index[s] for s in seq_dict[k]])
    keys_list.append(k)    
seq_msa = np.array(seq_msa)

## Separate the key lists of the training and testing dataset and save them.
training_keys_list = keys_list[0:33972]
testing_keys_list = keys_list[33972:]

# training keys
with open("./{}/keys_list.pkl".format(output_folder), 'wb') as file_handle:
    pickle.dump(training_keys_list, file_handle)

# testing keys
with open("./{}/t_keys_list.pkl".format(output_folder), 'wb') as file_handle:
    pickle.dump(testing_keys_list, file_handle)

## Separate seq_msa data.
training_seq_msa = seq_msa[0:33972,:]
testing_seq_msa = seq_msa[33972:,:]

## reweighting sequences
# note: only reweighting sequences that are used for training.
# seq_msa.shape[0]: number of sequences
# seq_msa.shape[1]: sequence length
seq_weight = np.zeros(training_seq_msa.shape)
for j in range(training_seq_msa.shape[1]):
    aa_type, aa_counts = np.unique(training_seq_msa[:,j], return_counts = True)
    num_type = len(aa_type)
    aa_dict = {}
    for a in aa_type:
        aa_dict[a] = aa_counts[list(aa_type).index(a)]
    for i in range(training_seq_msa.shape[0]):
        seq_weight[i,j] = (1.0/num_type) * (1.0/aa_dict[seq_msa[i,j]])
tot_weight = np.sum(seq_weight)
seq_weight = seq_weight.sum(1) / tot_weight 
with open("./{}/seq_weight.pkl".format(output_folder), 'wb') as file_handle:
    pickle.dump(seq_weight, file_handle)

# testing sequence weight
# All testing sequences have an equal weight.
t_seq_weight = np.ones(testing_seq_msa.shape[0]) / testing_seq_msa.shape[0]
t_seq_weight = t_seq_weight.astype(np.float32)
with open("./{}/t_seq_weight.pkl".format(output_folder), 'wb') as file_handle:
    pickle.dump(t_seq_weight, file_handle)

## remove positions where too many sequences have gaps
pos_idx = []
for i in range(training_seq_msa.shape[1]):
    if np.sum(seq_msa[:,i] == 0) <= seq_msa.shape[0]*0.2:
        pos_idx.append(i)
with open("./{}/seq_pos_idx.pkl".format(output_folder), 'wb') as file_handle:
    pickle.dump(pos_idx, file_handle)
    
training_seq_msa = training_seq_msa[:, np.array(pos_idx)]
testing_seq_msa = testing_seq_msa[:, np.array(pos_idx)]
with open("./{}/seq_msa.pkl".format(output_folder), 'wb') as file_handle:
    pickle.dump(training_seq_msa, file_handle)
with open("./{}/t_seq_msa.pkl".format(output_folder), 'wb') as file_handle:
    pickle.dump(testing_seq_msa, file_handle)

## change aa numbering into binary
K = 21 ## num of classes of aa
D = np.identity(K)
training_num_seq = training_seq_msa.shape[0]
testing_num_seq = testing_seq_msa.shape[0]
len_seq_msa = training_seq_msa.shape[1]
training_seq_msa_binary = np.zeros((training_num_seq, len_seq_msa, K))
testing_seq_msa_binary = np.zeros((testing_num_seq, len_seq_msa, K))
for i in range(training_num_seq):
    training_seq_msa_binary[i,:,:] = D[training_seq_msa[i]]
for j in range(testing_num_seq):
    testing_seq_msa_binary[j,:,:] = D[testing_seq_msa[j]]

## seq_msa_binary and t_seq_msa_binary are two major files that will be used for training and analysis.
# training binary
with open("./{}/seq_msa_binary.pkl".format(output_folder), 'wb') as file_handle:
    pickle.dump(training_seq_msa_binary, file_handle)

# testing binary
with open("./{}/t_seq_msa_binary.pkl".format(output_folder), 'wb') as file_handle:
    pickle.dump(testing_seq_msa_binary, file_handle)

