import pandas as pd
from Bio.Blast.Applications import NcbiblastxCommandline
import os
import io
import numpy as np
from Bio import Align
import argparse
import numpy as np
from Bio import SeqIO
from Bio import SearchIO
from Bio import Align
import re


def get_species_row(species_name, db_path_folder):
    """
    Reads the 'genomes_assemblies_CoreGenes_final.txt' file (without a header) into a DataFrame 
    and returns the row(s) where the first column matches the given species name.

    Parameters:
    - species_name (str): The species name to look for in the first column.

    Returns:
    - pd.DataFrame: A DataFrame containing the matching row(s) if found, or an empty DataFrame if not found.
    """
    file_path = db_path_folder+"genomes_assemblies_CoreGenes_final.txt"  # File path to the tab-separated file
    
    # Read the tab-separated file into a DataFrame without headers
    try:
        dataframe = pd.read_csv(file_path, sep="\t", header=None)
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        return pd.DataFrame()  # Return an empty DataFrame
    except pd.errors.EmptyDataError:
        print(f"Error: File '{file_path}' is empty.")
        return pd.DataFrame()  # Return an empty DataFrame
    except pd.errors.ParserError:
        print(f"Error: Unable to parse file '{file_path}'. Check the format.")
        return pd.DataFrame()  # Return an empty DataFrame

    # Check if the DataFrame is empty
    if dataframe.empty:
        print("Error: The DataFrame is empty. No data to search.")
        return dataframe

    # Check if the first column exists
    if dataframe.columns.size == 0:
        print("Error: The DataFrame has no columns. Cannot perform search.")
        return pd.DataFrame()  # Return an empty DataFrame
    
    # Rename columns assuming there are three columns in the file
    dataframe.columns = ["Species", "Reference_genome", "Core genes found in the reference genome"]

    # Filter the DataFrame for rows where the "Species" column matches the species name
    result = dataframe[dataframe["Species"].str.contains(species_name, case=False, na=False)]
    #print(result)
    result2 = result.iloc[0, :].to_frame().T
    return result2


# Function to parse a FASTA file and keep full headers with their identifiers
def extract_sequences_to_dict(fasta_file):
    sequences_dict = {}
    current_key = None  # To track the current header key

    with open(fasta_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):  # Header line
                current_key = line[1:].split()[0]  # Remove '>' and extract the first part of the header
                sequences_dict[current_key] = ""  # Initialize the key with an empty sequence
            elif current_key:  # Sequence lines
                sequences_dict[current_key] += line  # Append the sequence to the current key
    
    return sequences_dict


def run_blastx_for_references(query_file, reference_df, db_path_folder, num_threads=5):
    """
    Run BLASTX for a given query file against multiple reference databases specified 
    in the DataFrame's second column, and return the combined results as a DataFrame.

    Parameters:
    - query_file (str): Path to the query FASTA file.
    - reference_df (pd.DataFrame): DataFrame containing reference genome names in the second column.
    - num_threads (int): Number of threads to use for BLASTX. Defaults to 5.
    - evalue_threshold (float): E-value threshold for BLASTX. Defaults to 1e-10.
    - max_target_seqs (int): Maximum number of target sequences to report. Defaults to 1000.

    Returns:
    - pd.DataFrame: A DataFrame containing the combined BLASTX results for all reference genomes.
    """
    # Initialize an empty DataFrame to hold the results
    all_results_df = pd.DataFrame()
    
    # Iterate over the reference names in the second column of the DataFrame
    for index, row in reference_df.iterrows():
        reference_genome_name = row[1]  # Assuming the second column has the reference genome names
        
        # Define the database path (assuming the db name is the same as the reference genome name)
        db_path = db_path_folder+"blast_databases/"+reference_genome_name  # You may need to modify this if the path is different
        
        # Define the BLASTX command
        blastx_cline = NcbiblastxCommandline(
            query=query_file,              # Path to the query file
            db=db_path,                    # Path to the BLAST database
            outfmt="6 std qseq sseq slen", # Output format with query and subject sequences
            max_target_seqs=100,           # Maximum number of target sequences to report
            evalue=1e-10,                 # Reduced gap extend penalty
            num_threads=num_threads,       # Number of threads for parallel processing
            soft_masking=True,             # Enable soft masking for low-complexity regions
        )
        
        # Run the BLASTX command and capture the output
        print(f"Running BLASTX for reference: {reference_genome_name}")
        stdout, stderr = blastx_cline()
        
        if stderr:
            print(f"Error running BLASTX for {reference_genome_name}: {stderr}")
        else:
            # Read the BLASTX output into a DataFrame
            result_df = pd.read_csv(io.StringIO(stdout), sep="\t", header=None)
            
            # Add a column to indicate the reference genome
            result_df['Reference_Genome'] = reference_genome_name
            
            # Append the result to the combined DataFrame
            all_results_df = pd.concat([all_results_df, result_df], ignore_index=True)
    
    return all_results_df
    
# Function to read a FASTA file and create a dictionary with identifiers and gene names
def parse_fasta_to_dict(fasta_file):
    sequences_dict = {}
    # Regular expression to extract the identifier and gene name
    #header_pattern = re.compile(r'^>lcl\|([^\s]+).*\[gene=([^\]]+)\]')
    header_pattern = re.compile(r'^>lcl\|([^\s]+).+?\[gene=([^\]]+)\]')
    with open(fasta_file, 'r') as file:
        for line in file:
            if line.startswith('>'):  # Identify header lines
                match = header_pattern.search(line)
                if match:
                    #print(match.group)
                    identifier = match.group(1)  # Extract the unique identifier
                    gene_name = match.group(2)  # Extract the gene name
                    # Store in the dictionary
                    sequences_dict[identifier] = gene_name
    #print(sequences_dict)
    return sequences_dict
#Fuction for reading the fasta file 
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


# Function to get length from dictionary by matching name to a substring of the keys
def get_length(name, sequences_dict):
    # Iterate over the keys in the dictionary
    for key in sequences_dict.keys():
        if name in key:  # Check if the name is a substring of the key
            return sequences_dict.get(key, 0)  # Return the corresponding value if found
    return 0  # Default to 0 if no match is found


# Function to generate comments based on conditions in each row
def generate_comment(row):
    comments = []
    
    # Condition 1: If differences is [] and total_coverage <= 99
    if row['differences'] == [] and row['total_coverage'] <= 99:
        comments.append("Not the full protein")
    
    # Iterate over each element in the list "differences"
    for diff in row['differences']:
        for difference in diff: 
            if isinstance(difference, str) and "At position 1:" in difference:  # Check if the element contains "At position 1:"
                comments.append("Different first amino acid")
            # Check if the element matches "At position" with a number other than 1 and contains '->'
            elif isinstance(difference, str) and re.search(r"At position (\d+):.*->", difference):
                match = re.search(r"At position (\d+):", difference)
                if match:
                    position = int(match.group(1))  # Extract the position number
                    if position != 1:
                        comments.append("Amino acid mutation")
            
            elif "STOP CODON" in difference:
                comments.append("Early Stop Codon")
            elif re.search(r"deletion of", difference):
                comments.append("Deletion")
            elif re.search(r"insertion of", difference):
                comments.append("Insertion")  # Corrected this line by closing the parenthesis
    
    # Combine all generated comments with "; " as a separator, ensuring unique comments
    unique_comments = list(set(comments))  # Convert to set to remove duplicates and back to list
    return "; ".join(unique_comments)

# Function to find the matching gene name from the dictionary
def find_gene_name(identifier, genes_dict):
    # Iterate over the keys in the dictionary
    for key in genes_dict.keys():
        if key in identifier:  # Check if the key is a substring of the identifier
            return genes_dict[key]
    return 'Unknown'  # Default if no match is found

# Function to filter species rows based on exact match or partial match
def filter_species(row, species_short):
    # Normalize the species name by replacing spaces with underscores
    species_value = str(row['Species']).replace(' ', '_')  # Replace spaces with underscores
    # Check if the species is either exactly or partially matching the normalized Species column
    return species_short in species_value or species_value in species_short

# Function to extract positions and match mutations
def find_mutations(row, mutations_df, difference_number):
    positions = []
    # Ensure that 'differences' is a list of lists
    if isinstance(row['differences'], list):
        for diff in row['differences']:
            if isinstance(diff, list):
                for difference in diff:
                    if isinstance(difference, str):
                        match = re.search(r"At position (\d+):", difference)
                        if match:
                            positions.append(int(match.group(1)))
    # Find mutations within ±20 amino acid distance, ensuring positions are valid
    matched_mutations = []
    if positions:  # Proceed only if positions list is not empty
        for pos in positions:
            #print (pos)
            # Find all mutations that are within 20 amino acids above or below the position
            # and have the same Gene_ID as the input row
            nearby_mutations = mutations_df[
                (mutations_df['Gene_ID'] == row['Gene_Name']) &     # Check for matching Gene_ID
                mutations_df['Codon_pos'].notna() &               # Exclude NaN values
                (mutations_df['Codon_pos'] != 0) &  
                mutations_df['Codon_pos'].between(pos - difference_number, pos + difference_number)
            ]['Mutations_Resistance'].tolist()
            #print(nearby_mutations)
            matched_mutations.extend(nearby_mutations)
            #print(nearby_mutations)
    #return (matched_mutations)
    # Return the matched mutations as a semicolon-separated string, or indicate no matches found
    return "; ".join(matched_mutations) if matched_mutations else "No known mutations found"



# Define the main processing function
def process_blast_results(species, df, query_file, reference_df, difference_number, db_path_folder, difference_from_start = 20, identity_threshold=90, coverage_threshold=50):
    """
    Process BLAST results for a given species, query file, and reference data.
    
    Parameters:
    - species (str): The species name to analyze.
    - df (pd.DataFrame): The BLAST results DataFrame.
    - query_file (str): The path to the query FASTA file.
    - reference_df (pd.DataFrame): A one-row DataFrame with the reference genome name in the second column.
    - difference_number (int): A threshold number for finding mutations.

    Returns:
    - pd.DataFrame: A DataFrame containing the results of the processing.
    """
    # Extract the reference genome name from the reference DataFrame
    reference_genome_name = db_path_folder+"genomes_all/filtered_"+reference_df.iloc[0, 1]+"_translated_cds.faa"  # Assumes the second column has the reference genome name
    print(reference_genome_name)
    # Read sequences from FASTA files
    sequences_dict = read_fasta(query_file)
    genes_dict = parse_fasta_to_dict(reference_genome_name)
    reference_sequences = extract_sequences_to_dict(reference_genome_name)
    

    # Calculate additional columns
    df['contig_length'] = df[0].apply(lambda name: get_length(name, sequences_dict))
    df['total_coverage'] = (df[3] / df[14]) * 100
    df['partial_coverage'] = np.where(df[6] == 1, ((df[3] / (df[14] - df[3])) * 100), (df[3] / df[14]) * 100)
    df.loc[df['partial_coverage'] > 100, 'partial_coverage'] = 100
    # Filter the DataFrame
    
    df[[2, 6,7,8,10,11]] = df[[2, 6,7,8,10,11]].astype(float)
    df['partial_coverage'] = df['partial_coverage'].astype(float)
    #Change the subject query with the reference sequence
    df[13] = df[1].map(reference_sequences)
    #print(df)

    filtered_df = df[(df.iloc[:, 8] <= difference_from_start) | (df.iloc[:, 6] <= 10) | (df.iloc[:, 7] <= 10) | 
                     (df.iloc[:, 6] >= df.loc[:, "contig_length"]-10) | (df.iloc[:, 7] >= df.loc[:, "contig_length"]-10)]

    filtered_df2 = filtered_df[((filtered_df.iloc[:, 2] >= identity_threshold) & 
                                (filtered_df.loc[:, 'total_coverage'] >= coverage_threshold)) | 
                               ((filtered_df.iloc[:, 6] <= 10) & (filtered_df.iloc[:, 2] >= identity_threshold)) |
                               ((filtered_df.iloc[:, 7] <= 10) & (filtered_df.iloc[:, 2] >= identity_threshold)) |
                               ((filtered_df.iloc[:, 6] >= filtered_df.loc[:, "contig_length"]-10) & 
                                (filtered_df.iloc[:, 2] >= identity_threshold)) | 
                               ((filtered_df.iloc[:, 7] >= filtered_df.loc[:, "contig_length"]-10) & 
                                (filtered_df.iloc[:, 2] >= identity_threshold))]
    filtered_df2 = filtered_df2.reset_index(drop=True)

    # Identify non-functional proteins
    not_functional_prot = df[((df.iloc[:, 2] < identity_threshold) & (df.loc[:, 'total_coverage'] < coverage_threshold)) | 
                             ((df.iloc[:, 6] == 1) & (df.iloc[:, 2] >= identity_threshold)) |
                             ((df.iloc[:, 7] == 1) & (df.iloc[:, 2] >= identity_threshold)) |
                             ((df.iloc[:, 6] == df.loc[:, "contig_length"]) & (df.iloc[:, 2] >= identity_threshold)) | 
                             ((df.iloc[:, 7] == df.loc[:, "contig_length"]) & (df.iloc[:, 2] >= identity_threshold)) |
                             (df.iloc[:, 8] > difference_from_start)]

    # Sort and filter duplicates
    filtered_df3 = filtered_df2.sort_values(by=[10, 11, 13], ascending=[True, False, False])
    filtered_df4 = filtered_df3.drop_duplicates(subset=1, keep='first')

    # Add differences column
    filtered_df4['differences'] = 0
    filtered_df4['differences'] = filtered_df4['differences'].astype('object')
    filtered_df4['total_coverage'] = filtered_df4['total_coverage'].astype(float)

    # Alignment process and mutation identification
    for index, row in filtered_df4.iterrows():
        query_id = row[0]
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'  # Needleman-Wunsch for global alignment
        aligner.match_score = 1
        aligner.mismatch_score = -1
        aligner.open_gap_score = -2
        aligner.extend_gap_score = -0.5
        seq1= row[12]
        seq2= row[13]
        alignments = aligner.align(row[12], row[13])
        alignment = alignments[0]
        # Reconstruct the aligned query and target sequences with gaps
        aligned_query = []
        aligned_target = []

        # Extract the aligned blocks
        query_pos, target_pos = alignment.aligned

        # Keep track of positions in both sequences
        query_index = 0
        target_index = 0

        # Iterate over the aligned blocks
        for query_block, target_block in zip(query_pos, target_pos):
            # Handle any gaps before the current block
            while query_index < query_block[0] or target_index < target_block[0]:
                if query_index < query_block[0]:
                    aligned_query.append(seq1[query_index])
                    query_index += 1
                    aligned_target.append('-')  # Gap in target
                elif target_index < target_block[0]:
                    aligned_target.append(seq2[target_index])
                    target_index += 1
                    aligned_query.append('-')  # Gap in query

            # Handle the current aligned block
            for i in range(query_block[1] - query_block[0]):
                aligned_query.append(seq1[query_index])
                aligned_target.append(seq2[target_index])
                query_index += 1
                target_index += 1

        # Add trailing gaps after the final aligned block
        while query_index < len(seq1):
            aligned_query.append(seq1[query_index])
            query_index += 1
            aligned_target.append('-')

        while target_index < len(seq2):
            aligned_target.append(seq2[target_index])
            target_index += 1
            aligned_query.append('-')

        # Convert lists to strings
        aligned_query_str = ''.join(aligned_query)
        aligned_target_str = ''.join(aligned_target)

        query_aligned = aligned_query_str
        subject_aligned = aligned_target_str
        #print(query_aligned)
        # Identify differences
        differences = []
        prev_event = None
        prev_position = None
        start_position = None
        continuous_change = ""
        change = ""  # Initialize change variable with an empty string

        for i, (q, s) in enumerate(zip(query_aligned, subject_aligned)):
            if q != s:
                # Determine the type of event
                if (q == "-") and (s != "-"):
                    event = "deletion"
                    change = s
                elif (q != "-") and (s == "-"):
                    event = "insertion"
                    change = q
                else:
                    event = "substitution"
                    change = f"{s} -> {q}"

                # Check if this event is continuous
                if event == prev_event and i == prev_position + 1 and event != "substitution":
                    # Extend the range for a continuous event (not substitutions)
                    continuous_change += change
                else:
                    # Finalize the previous event, if any
                    if prev_event is not None:
                        if prev_event != "substitution":
                            # Report the continuous range
                            end_position = prev_position + 1
                            differences[-1] = f"At positions {start_position}-{end_position}: {prev_event} of {continuous_change}"
                        else:
                            # Report substitution
                            differences[-1] = f"At position {prev_position + 1}: {differences[-1].split(': ')[-1]}"
                    
                    # Start a new event
                    start_position = i + 1
                    if event != "substitution":
                        continuous_change = change
                        differences.append(f"At position {i + 1}: {event}")
                    else:
                        differences.append(f"At position {i + 1}: {event} of {change}")

                prev_event = event
                prev_position = i
            else:
                # Finalize any ongoing continuous event when alignment matches
                if prev_event is not None and prev_event != "substitution":
                    end_position = prev_position + 1
                    differences[-1] = f"At positions {start_position}-{end_position}: {prev_event} of {continuous_change}"
                prev_event = None
                continuous_change = ""

        # Finalize the last event, if necessary
        if prev_event is not None and prev_event != "substitution":
            end_position = prev_position + 1
            differences[-1] = f"At positions {start_position}-{end_position}: {prev_event} of {continuous_change}"
        elif prev_event == "substitution":
            differences[-1] = f"At position {prev_position + 1}: {differences[-1].split(': ')[-1]}"

        value_for_new_column = differences if differences else []
        wrapped_lists = [[lst] for lst in value_for_new_column]
        filtered_df4.at[index, 'differences'] = wrapped_lists

    # Add Gene_Name and identify missing genes
    filtered_df4['Gene_Name'] = filtered_df4[1].apply(lambda name: find_gene_name(name, genes_dict))
    all_genes = set(genes_dict.values())
    genes_in_df = set(filtered_df4['Gene_Name'])
    missing_genes = all_genes - genes_in_df

    # Create DataFrame for missing genes
    missing_rows = pd.DataFrame({
        0: [f'Missing_{i}' for i in range(1, len(missing_genes) + 1)],  # Fill identifiers as needed
        1: [None] * len(missing_genes),  # Placeholder for missing values in column 1
        2: [None] * len(missing_genes),
        3: [None] * len(missing_genes),
        4: [None] * len(missing_genes),
        5: [None] * len(missing_genes),
        6: [None] * len(missing_genes),
        7: [None] * len(missing_genes),
        8: [None] * len(missing_genes),
        9: [None] * len(missing_genes),
        'contig_length': [0] * len(missing_genes),  # Default values for missing rows
        'total_coverage': [0.0] * len(missing_genes),  # Example default values
        'partial_coverage': [0.0] * len(missing_genes),  # Default partial coverage
        'differences': [[] for _ in range(len(missing_genes))],  # Empty differences
        'Gene_Name': list(missing_genes),
        'Comments': ['Not functional protein / Missing protein' for _ in range(len(missing_genes))]
    })

    # Filter rows where "differences" is [] and "total_coverage" is >= 90
    same_proteins = filtered_df4[(filtered_df4['differences'].apply(lambda x: x == [])) & 
                                 (filtered_df4['total_coverage'] >= 90) & (filtered_df4[2] == 100.0)]
    different_proteins = filtered_df4[~((filtered_df4['differences'].apply(lambda x: x == [])) & 
                                       (filtered_df4['total_coverage'] >= 90))]
    different_proteins['Comments'] = different_proteins.apply(generate_comment, axis=1)

    # Append missing rows to the original DataFrame
    missing_rows = missing_rows.astype(str)
    different_proteins = pd.concat([different_proteins, missing_rows], ignore_index=True)
    #different_proteins = pd.concat([different_proteins, missing_rows], ignore_index=True)
	
    new_names = {
        0: 'Contig ID',
        1: 'ID',
        2: "Identity",
        3: "Length Alignment",
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
        14: "Length of Reference Protein", 
        "Reference_Genome" : "Reference_Genome",
        'contig_length': 'contig_length',
        'total_coverage': 'total_coverage',
        'partial_coverage': 'partial_coverage',
        'differences': 'differences',
        'Gene_Name': 'Gene_Name',
        'Comments' : 'Comments'
    }
    different_proteins = different_proteins.rename(columns=new_names)
    
    
    # Process mutations details
    mutations_details = pd.read_csv(db_path_folder+'all_detailed/full_list_mutations.csv', sep=',')
    mutations_details['Codon_pos'] = pd.to_numeric(mutations_details['Codon_pos'], errors='coerce')
    mutations_details['Codon_pos'] = mutations_details['Codon_pos'].fillna(0).astype(int)
    
    # Filter out 'Mycobacterium' species
    species_short = '_'.join(species.split(' ')[:2])
    mutations_details = mutations_details[~mutations_details['Species'].str.contains('Mycobacterium', na=False, case=False)]
    mutations_details['Normalized_Species'] = mutations_details['Species'].str.replace(' ', '_')

    # Create a mask for rows that match the species criteria
    match_mask = mutations_details['Normalized_Species'].notna() & mutations_details.apply(
        lambda row: filter_species(row, species_short), axis=1)
    mutations_details_Species = mutations_details[match_mask]
    mutations_details_OtherSpecies = mutations_details[~match_mask]
	
	
    # Step 1: Extract the first two words of the species name and format with an underscore
    formatted_species = "_".join(species.split()[:2])

    # Step 2: Dynamically create the column name based on the formatted species
    column_name = f"{formatted_species} Known Mutations                                                                                                                                                                                                                                                                               ."


    # Apply the function to create the 'Species_Mutations' column
    different_proteins[column_name] = different_proteins.apply(
        lambda row: find_mutations(row, mutations_details_Species,difference_number), axis=1
    )

# Apply the function to create the 'Species_Mutations' column
    different_proteins['Other_Species_Known_Mutations']=different_proteins.apply(
        lambda row: find_mutations(row, mutations_details_OtherSpecies,difference_number), axis=1
    )

    different_proteins2 = different_proteins.loc[:,['Contig ID', 'Gene_Name', "ID","Comments",
	"differences",column_name,"Other_Species_Known_Mutations",
	"Identity", "partial_coverage", "Length of Reference Protein",
	"Start", "Stop", "Start In Reference Protein", "Stop In Reference Protein", 
	"Length Alignment", "Mismatches", "Gaps", "E-value","BitScore","Reference_Genome",
	"Genome Protein Sequence", "Protein Reference Sequence"]] 
    different_proteins2.rename(columns={'ID': 'Core Gene ID', "Gene_Name" : "Gene Name",
    			"Start" : "Query Start", "Stop" : "Query Stop",
    			"differences" : "Differences from Reference Protein",
    			"Start In Reference Protein" : "Reference Start", 
    			"Stop In Reference Protein" : "Reference Stop",
    			"Other_Species_Known_Mutations" : "Other Species Known Mutations                                                                                                                                                                                                                                                                               .",
    			"Identity" : "Identity", 
    			"partial_coverage": "Coverage",
    			"Length Alignment" : "Alignment Length",
    			"Genome Protein Sequence" : "Genome Protein Sequence                                                                                                                                                                                                                                                                                                                                                                                                         .",
    			"Protein Reference Sequence" : "Protein Reference Sequence                                                                                                                                                                                                                                                                                                                                                                                                         ."}, inplace=True) 
    				# Specify the column you want to keep in its original data type
    column_to_exclude = ["Identity","Coverage",
    				 "E-value","BitScore"]

    different_proteins2['Core Gene ID'] = different_proteins2['Core Gene ID'].str.replace(r'^lcl\|', '', regex=True)
    different_proteins2['Core Gene ID'] = different_proteins2['Core Gene ID'].str.replace(r'_\d+$', '', regex=True)   # Remove '_123' at the end
    # Convert all columns to strings without decimals, except the specified columns
    different_proteins2 = different_proteins2.apply(lambda x: x.fillna(0).astype(int).astype(str) if x.name not in column_to_exclude and x.dtype != 'object' else x)

    return different_proteins2



