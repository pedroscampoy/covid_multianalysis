import os
import sys
import re
import argparse
import subprocess
import pandas as pd
import numpy as np
import shutil


def get_arguments():

    parser = argparse.ArgumentParser(prog = 'snp_covid19.py', description = 'Pipeline to do variant calling with Sars-CoV-2 samples')

    input_group = parser.add_argument_group('Input', 'Input parameters')

    input_group.add_argument('-i', '--input', dest="input_dir", type=str, required=True, help="Required. Input dir with the fastq files")

    output_group = parser.add_argument_group('Output', 'Output parameters')

    output_group.add_argument('-o', '--output', type=str, required=True, help="Required.Output directory to extract all results")

    arguments = parser.parse_args()

    return arguments

args = get_arguments()

def check_create_dir(path):
    #exists = os.path.isfile(path)
    #exists = os.path.isdir(path)
    if os.path.exists(path):
        pass
    else:
        os.mkdir(path)

out_compare_dir = os.path.join(args.output, "Compare")
input_dir = os.path.abspath(args.input_dir)

check_create_dir(out_compare_dir)

# 2 parts: presence matrix, pairwise

#1 - Presence matrix

run_id = input_dir.split("/")[-1]

def blank_database():
    new_pandas_ddtb = pd.DataFrame(columns=['Position','N', 'Samples'])
    return new_pandas_ddtb

final_ddbb = blank_database()
sample_filter_list = []
all_samples = 0
new_samples = 0
for r,d,f in os.walk(input_dir):
    for file in f:
        if file.endswith(".tsv"):
            all_samples = all_samples + 1
            positions_shared = []
            positions_added = []
            file_name = file.split(".")[0]
            print(file_name)
            file_path = os.path.join(r, file)
            print(file_path)
            read_file = pd.read_csv(file_path, sep="\t") #read the file
            file_pass = read_file[read_file["PASS"] == True] ###
            final_file = file_pass[~file_pass.ALT.str.startswith("-")&~file_pass.ALT.str.startswith("+")] ###
            for position in final_file["POS"].unique(): #access to the snps of each sample
                if position not in final_ddbb["Position"].values:
                    positions_added.append(int(position))
                    
                    new_row = len(final_ddbb.index)
                    final_ddbb.loc[new_row, "Position"] = int(position)
                    final_ddbb.loc[new_row, "Samples"] = file_name
                    final_ddbb.loc[new_row, "N"] = int(1)
                    final_ddbb.loc[new_row, file_name] = str(1)
                else:
                    positions_shared.append(int(position))
                    
                    index_position = final_ddbb.index[final_ddbb["Position"] == int(position)][0]
                    
                    number_samples_with_position = final_ddbb.loc[index_position, "N"]
                    names_samples_with_position = final_ddbb.loc[index_position, "Samples"]
                    new_names_samples = names_samples_with_position + "," + file_name
                    
                    final_ddbb.loc[index_position, "N"] = number_samples_with_position + 1
                    final_ddbb.loc[index_position, "Samples"] = new_names_samples
                    final_ddbb.loc[index_position, file_name] = str(1)
                
final_ddbb = final_ddbb.fillna(0)      
final_ddbb["Position"] = final_ddbb["Position"].astype(int)
final_ddbb["N"] = final_ddbb["N"].astype(int)

file_presence = run_id + ".presence.tsv"
file_presence_def = os.path.join(out_compare_dir, file_presence)

final_ddbb.to_csv(file_presence_def, sep="\t", index=False)

#2 - Pairwise

from sklearn.metrics import pairwise_distances, accuracy_score

def compare_snp_columns(sample1, sample2, df):
    jaccard_similarity = accuracy_score(df[sample1], df[sample2]) #similarities between colums
    hamming_similarity = 1 - jaccard_similarity #disagreements between colums
    snp_distance = int(hamming_similarity * (len(df.index)+1))
    return snp_distance

def snp_distance_pairwise(dataframe, output_file):
    with open(output_file, "a") as f:
        for sample1 in dataframe.iloc[:,3:].columns:
            for sample2 in dataframe.iloc[:,3:].columns:
                if sample1 != sample2:
                    snp_distance = compare_snp_columns(sample1, sample2, dataframe)
                    line_distance = "%s\t%s\t%s\n" % (sample1, sample2, snp_distance)
                    f.write(line_distance)

def import_to_pandas(file_table, header=False, sep='\t'):
    if header == False:
        #exclude first line, exclusive for vcf outputted by PipelineTB
        dataframe = pd.read_csv(file_table, sep=sep, skiprows=[0], header=None)
    else:
        #Use first line as header
        dataframe = pd.read_csv(file_table, sep=sep, header=0)
    
    return dataframe

file_presence_again = run_id + ".presence.tsv"
file_presence_again_def = os.path.join(out_compare_dir, file_presence_again) # instead of using previous one, we upload again, as in jupyter there was an error if we don't upload de output file

presence_ddbb = import_to_pandas(file_presence_again_def, header=True)

file_pairwise = run_id + ".snps.pairwise.tsv"
file_pairwise_path = os.path.join(out_compare_dir, file_pairwise)

snp_distance_pairwise(presence_ddbb, file_pairwise_path)