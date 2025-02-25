# find_amr_genes_module.py

import pandas as pd
import numpy as np
import os
from Bio.Blast.Applications import NcbiblastxCommandline
import io
import subprocess

ratio_threshold = 0.2

# Get the directory where find_amr_genes_module.py is located
MODULE_DIR = os.path.dirname(os.path.abspath(__file__))

def get_blast_version():
    try:
        # Run the 'blastx -version' command and capture the output
        result = subprocess.run(["blastx", "-version"], capture_output=True, text=True)
        
        # Print the captured output
        if result.returncode == 0:
            return result.stdout.strip()
        else:
            return f"Error: {result.stderr.strip()}"
    except FileNotFoundError:
        return "BLASTX not found. Please ensure BLAST+ is installed and accessible in your PATH."


def run_blastx_to_dataframe(query_file, db_path, num_threads=5):
    """
    Run BLASTX with the specified parameters and return the results as a DataFrame.

    Parameters:
    - query_file: Path to the query FASTA file.
    - db_path: Path to the BLAST database.
    - output_file: Path to the output file.
    - num_threads: Number of threads to use for BLASTX.

    Returns:
    - A pandas DataFrame with BLASTX results.
    """
    blastx_cline = NcbiblastxCommandline(
        query=query_file,
        db=db_path,
        outfmt="6 std qseq sseq slen",
        num_alignments =10000,
        evalue=1e-10,  # Set the e-value threshold
        num_threads=num_threads,
        soft_masking=True, 
    )
    
    stdout, stderr = blastx_cline()
    
    # Check for errors in stderr
    if stderr:
        raise Exception(f"Error running BLASTX: {stderr}")

    # Read the BLASTX results into a DataFrame directly from stdout
 #   column_names = [
 #       "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
 #       "qstart", "qend", "sstart", "send", "evalue", "bitscore", 
 #       "qseq", "sseq", "slen"
 #   ]

    # Use io.StringIO to read the stdout string into a DataFrame
    blast_df = pd.read_csv(io.StringIO(stdout), sep="\t", header=None)
    
    return blast_df


def save_dataframe_to_csv(dataframe, output_csv_path):
    """
    Save the given DataFrame to a CSV file.

    Parameters:
    - dataframe: The pandas DataFrame to save.
    - output_csv_path: Path to the output CSV file.
    """
    dataframe.to_csv(output_csv_path, index=False)


def extract_sequences_to_dict(fasta_file):
    """
    Parse a FASTA file and keep full headers with their identifiers.
    Returns a dictionary mapping sequence IDs to their full sequences.
    """
    sequences_dict = {}
    current_key = None

    with open(fasta_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):  # Header line
                current_key = line[1:].split()[0]  # Remove '>' and extract the first part of the header
                sequences_dict[current_key] = ""  # Initialize the key with an empty sequence
            elif current_key:  # Sequence lines
                sequences_dict[current_key] += line  # Append the sequence to the current key
    
    return sequences_dict


def process_blast_results(protein_annotation_file, assembly_blast_results, file_path, protein_start_filter=50, identity_threshold=70, coverage_threshold=70):
    """
    Process BLAST results, filter based on thresholds, and annotate.

    Parameters:
    - protein_annotation_file: Path to the protein annotation CSV file.
    - file_path: Path to the FASTA file for length calculations.
    - file_output: Path to the output CSV file.
    - protein_start_filter: Threshold for protein start position filter.
    - identity_threshold: Threshold for sequence identity.
    - coverage_threshold: Threshold for coverage percentage.

    Returns:
    - None
    """
    protein_annotation = pd.read_csv(protein_annotation_file, delimiter='\t')
    df = assembly_blast_results

    # Extract reference sequences from the AMR proteins FASTA file
    amr_proteins_fasta = os.path.join(MODULE_DIR, "databases", "all", "nr100_all_proteins_amr.fasta")
    reference_sequences = extract_sequences_to_dict(amr_proteins_fasta)
    
    # Map the reference sequences to column 13
    df[13] = df[1].map(reference_sequences)

    def read_fasta(file_path):
        sequences = {}
        with open(file_path, 'r') as file:
            entry_name = None
            sequence_length = 0
            for line in file:
                line = line.strip()
                if line.startswith('>'):
                    if entry_name is not None:
                        sequences[entry_name] = sequence_length
                    entry_name = line[1:]
                    sequence_length = 0
                else:
                    sequence_length += len(line)
            if entry_name is not None:
                sequences[entry_name] = sequence_length
        return sequences

    sequences_dict = read_fasta(file_path)
    df['total_coverage'] = (df[3] / df[14]) * 100

    def get_length(name):
        for key in sequences_dict.keys():
            if name in key:
                return sequences_dict.get(key, 0)
        return 0

    # Create the new column using numpy.where()
    df['partial_coverage'] = np.where((df[6] == 1) | (df[7] == 1), 
                                  (df[3] / (df[14] - df[3])) * 100, 
                                  (df[3] / df[14]) * 100)

    df.loc[df['partial_coverage'] > 100, 'partial_coverage'] = 100
    df['contig_length'] = df[0].apply(get_length)

    df[[2, 6,7,8,10,11]] = df[[2, 6,7,8,10,11]].astype(float)
    df['partial_coverage'] = df['partial_coverage'].astype(float)
    
    filtered_df = df[(df.iloc[:, 8] <= protein_start_filter) | 
                     (df.iloc[:, 6] <= 10) | 
                     (df.iloc[:, 7] <= 10) | 
                     (df.iloc[:, 6] >= df.loc[:, "contig_length"]-10) | 
                     (df.iloc[:, 7] >= df.loc[:, "contig_length"]-10)]

    filtered_df2 = filtered_df[((filtered_df.iloc[:, 2] >= identity_threshold) & 
                                (filtered_df.loc[:, 'partial_coverage'] >= coverage_threshold)) | 
                               ((filtered_df.iloc[:, 6] <= 10) &   
                                (filtered_df.iloc[:, 2] >= identity_threshold)) |
                               ((filtered_df.iloc[:, 7] <= 10) &   
                                (filtered_df.iloc[:, 2] >= identity_threshold)) |
                               ((filtered_df.iloc[:, 6] >= filtered_df.loc[:, "contig_length"]-10) & 
                                (filtered_df.iloc[:, 2] >= identity_threshold)) | 
                               ((filtered_df.iloc[:, 7] >= filtered_df.loc[:, "contig_length"]-10) & 
                                (filtered_df.iloc[:, 2] >= identity_threshold))]
    filtered_df2 = filtered_df2.reset_index(drop=True)

    filtered_df2[[2, 10,11]] = filtered_df2[[2, 10,11]].astype(float)
    filtered_df2['partial_coverage'] = filtered_df2['partial_coverage'].astype(float)

    df1 = pd.DataFrame()
    for index, row in filtered_df2.iterrows():
        inner_break = False
        node = row.iloc[0]
        min_value = min(row.iloc[7], row.iloc[6])
        max_value = max(row.iloc[7], row.iloc[6])
        for index_2, row_2 in filtered_df2[filtered_df2[0] == node].iterrows():
            if row.equals(row_2):
                continue
            else:
                min_value_2 = min(row_2.iloc[7], row_2.iloc[6])
                max_value_2 = max(row_2.iloc[7], row_2.iloc[6])
                if (min_value >= max_value_2) or (max_value <= min_value_2):
                    continue
                elif (min_value >= min_value_2) & (min_value <= max_value_2) & (max_value >= max_value_2):
                    ratio = float((max_value_2 - min_value) / (max_value_2 - min_value_2))
                    if ratio >= ratio_threshold:
                        if row.loc["partial_coverage"] > row_2.loc["partial_coverage"] and row.iloc[2] > row_2.iloc[2]:
                            continue
                        elif row.loc["partial_coverage"] >= row_2.loc["partial_coverage"] and row.iloc[2] > row_2.iloc[2]:
                            continue
                        elif row.loc["partial_coverage"] > row_2.loc["partial_coverage"] and row.iloc[2] >= row_2.iloc[2]:
                            continue
                        # Case 1: First value is strictly greater AND second value is less than or equal
                        elif (row.iloc[2] > row_2.iloc[2] and row.iloc[10] <= row_2.iloc[10]) or (row.iloc[2] >= row_2.iloc[2] and row.iloc[10] < row_2.iloc[10]):
                            continue
                        # Case 1: First value is strictly greater AND second value is greater than or equal
                        elif (row.iloc[2] > row_2.iloc[2] and row.iloc[11] >= row_2.iloc[11]) or  (row.iloc[2] >= row_2.iloc[2] and row.iloc[11] > row_2.iloc[11]):
                            continue
                        elif (row.iloc[2] >= row_2.iloc[2] and row.iloc[11] >= row_2.iloc[11]-2):
                            continue
                        else:
                            inner_break = True
                            break
                    else:
                        continue
                elif ((min_value < min_value_2) & (max_value <= max_value_2) & (max_value >= min_value_2)) or \
            ((min_value <= min_value_2) & (max_value < max_value_2) & (max_value >= min_value_2)):
                    ratio = float((max_value - min_value_2) / (max_value_2 - min_value_2))
                    if ratio >= ratio_threshold:
                        if row.loc["partial_coverage"] > row_2.loc["partial_coverage"] and row.iloc[2] > row_2.iloc[2]:
                            continue
                        elif row.loc["partial_coverage"] >= row_2.loc["partial_coverage"] and row.iloc[2] > row_2.iloc[2]:
                            continue
                        elif row.loc["partial_coverage"] > row_2.loc["partial_coverage"] and row.iloc[2] >= row_2.iloc[2]:
                            continue
                        # Case 1: First value is strictly greater AND second value is less than or equal
                        elif (row.iloc[2] > row_2.iloc[2] and row.iloc[10] <= row_2.iloc[10]) or  (row.iloc[2] >= row_2.iloc[2] and row.iloc[10] < row_2.iloc[10]):
                            continue
                        # Case 1: First value is strictly greater AND second value is greater than or equal
                        elif (row.iloc[2] > row_2.iloc[2] and row.iloc[11] >= row_2.iloc[11]) or (row.iloc[2] >= row_2.iloc[2] and row.iloc[11] > row_2.iloc[11]):
                            continue
                        # Case 1: First value is strictly greater AND second value is less than or equal
                        elif (row.iloc[2] > row_2.iloc[2] and row.iloc[11] >= row_2.iloc[11]-2) or (row.iloc[2] >= row_2.iloc[2] and row.iloc[11] > row_2.iloc[11]-2):
                            continue
                        else:
                            inner_break = True
                            break
                elif (min_value>=min_value_2) & (min_value<=max_value_2) & (max_value<=max_value_2):
                	ratio=float((max_value-min_value)/(max_value_2-min_value_2))
                	if ((ratio>=ratio_threshold) and (row.loc["partial_coverage"]>=80)):
                		if (row.loc["partial_coverage"]>=row_2.loc["partial_coverage"]) and (row.iloc[2]>=row_2.iloc[2]):
                     #   print("Yeah")
                			continue
                    	
                		else:
                     #   print("No")
                			inner_break = True
                			break  # Break out of the inner loop
                	else:
                		continue		
        if not inner_break:
            row_df = pd.DataFrame([row])
            df1 = pd.concat([df1, row_df], ignore_index=True)

    # Sort by column 2 in descending order to prioritize the highest values
    df1 = df1.sort_values(by=2, ascending=False)

    # Drop duplicates based on the column you want to deduplicate (e.g., column 12), keeping the first (highest) occurrence
    df1 = df1.drop_duplicates(subset=[0,12], keep='first')

    df1['Comments'] = 'Complete gene based on threshold'
    df1.loc[((df1[6] == 1) | (df1[7] == 1)), 'Comments'] = 'Partial gene in the start of Contig'
    df1.loc[((df1[6] == df1['contig_length']) | (df1[7] == df1['contig_length'])), 'Comments'] = 'Partial gene in the end of Contig'
    df1['GeneImportance'] = "Major"
    
    new_names = {
        0: 'Contig ID',
        1: 'ID',
        2: "% Identity with reference protein",
        3: "Alignment Length",
        4: "Mismatches",
        5: "Gaps",
        6: "Start",
        7: "Stop",
        8: "Start In Reference Protein",
        9: "Stop In Reference Protein",
        10: "E-value",
        11: "BitScore",
        12: "Genome Protein Sequence",
        13: "Protein Reference Sequence",
        14: "Length of reference protein", 
        "partial_coverage": "% Coverage of reference protein",
        "contig_length": "Length of the Identified Contig"
    }
    df1 = df1.rename(columns=new_names)

    protein_annotation['ID'] = protein_annotation['ID'].astype(str)

    for index, row in df1.iterrows():
        value = row['ID']
        row_extra = protein_annotation[protein_annotation['ID'].apply(lambda x: x in value)]
        
        if not row_extra.empty:
            for column in row_extra.columns:
                new_column_name = f"{column}"
                df1.loc[index, new_column_name] = row_extra[column].values[0]
    
    df1.loc[df1.loc[:, "Database"] == "Card Database", 'GeneImportance'] = "Minor"
    #print(df1.loc[:, "Database"])
    df3 = df1.sort_values(by='GeneImportance')
    
    df4 = df3.loc[:,['Contig ID', 'gene_family', "product_name",'ID', "% Identity with reference protein",
               "% Coverage of reference protein",
              "Mismatches", "Gaps", "class", "subclass", "Comments",
               "Resistance Mechanism",
               "Start", "Stop", "Start In Reference Protein", "Stop In Reference Protein", 
               "genbank_nucleotide_accession", 
              "pubmed_reference", "Database", "Notes", "Alignment Length", 
            "Length of reference protein", "E-value","BitScore",
             "Genome Protein Sequence", "Protein Reference Sequence"]]  
    df4.rename(columns={'gene_family': 'Resistance Gene', 'ID': 'Gene ID', 
    			"product_name" : "Name of the product of the Resistance Gene",
    			"% Identity with reference protein" : "Identity",
               "% Coverage of reference protein" : "Coverage",
               "Comments" : "Comments based on the results of Blastx", 
    			"class": "Antibiotic Class", "subclass": "Antibiotic Sublasses from the Databases", 
    			"Start" : "Query Start", "Stop" : "Query Stop", 
    			"Start In Reference Protein" : "Reference Start", 
    			"Stop In Reference Protein" : "Reference Stop", 
    			"Notes" : "Notes for the gene in the Database",
    			"genbank_nucleotide_accession" : "Genbank Accession",  
    			"Length of reference protein" : "Reference length",
    			"pubmed_reference" : "Pubmed",
    			"Genome Protein Sequence" : "Genome Protein Sequence                                                                                                                                                                                                                                                                                                                                                                                                         .",
    			"Protein Reference Sequence" : "Protein Reference Sequence                                                                                                                                                                                                                                                                                                                                                                                                         ."}, inplace=True) 
    			# Specify the column you want to keep in its original data type
    column_to_exclude = ["Identity","Coverage",
    				 "E-value","BitScore"]


    # Convert all columns to strings without decimals, except the specified columns
    df4 = df4.apply(lambda x: x.astype(int).astype(str) if x.name not in column_to_exclude and x.dtype != 'object' else x)


    #df3.to_csv(file_output, index=False)
    #print(f"Final results have been saved to {file_output}.")
    return df4
