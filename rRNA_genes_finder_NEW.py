import pandas as pd
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
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
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Seq import Seq

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


def get_species_row(species_name, db_path):
    """
    Reads the 'genomes_assemblies_CoreGenes.txt' file (without a header) into a DataFrame 
    and returns the row(s) where the first column matches the given species name.

    Parameters:
    - species_name (str): The species name to look for in the first column.

    Returns:
    - pd.DataFrame: A DataFrame containing the matching row(s) if found, or an empty DataFrame if not found.
    """
    file_path = db_path+"genomes_assemblies_rRNA_final.txt"  # File path to the tab-separated file
    
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
    dataframe.columns = ["Species", "Reference_genome", "rRNA"]

    # Filter the DataFrame for rows where the "Species" column matches the species name
    result = dataframe[dataframe["Species"].str.contains(species_name, case=False, na=False)]

    # Remove the unwanted suffix from the "Reference_genome" column
    if not result.empty:
        result["Reference_genome"] = result["Reference_genome"].str.replace("_rna_from_genomic.fna", "", regex=False)
    # Function to extract the counts and strings inside parentheses for 16S, 23S, and 5S
    def extract_rRNA_info(row, rRNA_type):
    # Modify the regex to be non-greedy by using .*? to capture all matches correctly
        pattern = rf"{rRNA_type}.*?\(([^)]+)\)"
        matches = re.findall(pattern, row)
    
        if matches:
            return f"{len(matches)} ({', '.join(matches)})"  # Return count and strings in parentheses
        else:
            return '0 (None)'  # Return '0' if no match is found

    # Create new columns for 16S rRNA, 23S rRNA, and 5S rRNA
    result['16S rRNA'] = result['rRNA'].apply(lambda x: extract_rRNA_info(x, '16S'))
    result['23S rRNA'] = result['rRNA'].apply(lambda x: extract_rRNA_info(x, '23S'))
    result['5S rRNA'] = result['rRNA'].apply(lambda x: extract_rRNA_info(x, '5S'))
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
                current_key = line[5:].split()[0]  # Remove '>' and extract the first part of the header
                sequences_dict[current_key] = ""  # Initialize the key with an empty sequence
            elif current_key:  # Sequence lines
                sequences_dict[current_key] += line  # Append the sequence to the current key
    
    return sequences_dict


def run_blastn_for_references(query_file, reference_df, db_path, num_threads=5):
    """
    Run BLASTN for a given query file against multiple reference databases specified 
    in the DataFrame's second column, and return the combined results as a DataFrame.

    Parameters:
    - query_file (str): Path to the query FASTA file.
    - reference_df (pd.DataFrame): DataFrame containing reference genome names in the second column.
    - num_threads (int): Number of threads to use for BLASTN. Defaults to 5.
    - evalue_threshold (float): E-value threshold for BLASTN. Defaults to 1e-10.
    - max_target_seqs (int): Maximum number of target sequences to report. Defaults to 1000.

    Returns:
    - pd.DataFrame: A DataFrame containing the combined BLASTN results for all reference genomes.
    """
    # Assuming df4.iloc[:, 2] contains the single string as a list
    rna_list_string = reference_df.iloc[:, 2].tolist()[0]

    # Split the string into a list of individual items
    rna_list = rna_list_string.split(', ')

    rna_dict = {}

    for item in rna_list:
            # Split the string to separate the part outside and inside the parentheses
            value, full_id = item.split('(')

            # Remove the closing parenthesis from the full_id and strip any extra whitespace
            full_id = full_id.rstrip(')')

            if " |" in full_id:
                value1, value2 = full_id.split(" |")
            else:
                value1 = full_id
                value2 = None  # or any other default value if no separator is present

            # Add to the dictionary
            rna_dict[value1] = item
    # Initialize an empty DataFrame to hold the results
    all_results_df = pd.DataFrame()
    #print(all_results_df)
    # Iterate over the reference names in the second column of the DataFrame
    for index, row in reference_df.iterrows():
        reference_genome_name = row[1]  # Assuming the second column has the reference genome names
        genome_ref_path = db_path+"genomes/filtered_" + reference_genome_name + "_rna_from_genomic.fna"  # You may need to modify this if the path is different
       	#print(genome_ref_path)
        # Define the database path (assuming the db name is the same as the reference genome name)
        db_path = db_path+"blast_databases/" + reference_genome_name  # You may need to modify this if the path is different
         

        # Define the BLASTN command
        blastn_cline = NcbiblastnCommandline(
            query=query_file,
            db=db_path,
            outfmt="6 std qseq sseq slen",
            max_target_seqs=100,
            evalue=1e-10,
            num_threads=num_threads
        )
        
        # Run the BLASTN command and capture the output
        print(f"Running BLASTN for reference: {reference_genome_name}")
        stdout, stderr = blastn_cline()
        
        if stderr:
            print(f"Error running BLASTN for {reference_genome_name}: {stderr}")
        else:
            # Read the BLASTN output into a DataFrame
            result_df = pd.read_csv(io.StringIO(stdout), sep="\t", header=None)
            
            # Add a column to indicate the reference genome
            result_df['Reference_Genome'] = reference_genome_name
            
            # Append the result to the combined DataFrame
            all_results_df = pd.concat([all_results_df, result_df], ignore_index=True)
    
    # Step 1: Create a new column 'ID' by removing 'lcl|' from the second column (index 1)
    all_results_df.iloc[:, 1] = all_results_df.iloc[:, 1].str.replace(r'lcl|', '')
    
    reference_sequences = extract_sequences_to_dict(genome_ref_path)
    print(rna_dict)
    # Step 2: Create a new column 'rRNA_gene' based on the 'ID' column using `rna_dict`
    all_results_df['rRNA_gene'] = all_results_df.iloc[:, 1].map(rna_dict)
    
    all_results_df[13] = all_results_df[1].map(reference_sequences)

    print("all_resuls")
    print(all_results_df.iloc[:, 1])
    print(all_results_df)
    return all_results_df
    

def find_rRNA_differences(species,df4, blast_results, query_file, db_path, difference_number = 5, num_threads=5):
	sequences_dict = read_fasta(query_file)
	df = blast_results
	df['contig_length'] = df[0].apply(lambda name: get_length(name, sequences_dict))
	df['total_coverage'] = (df[3] / df[14]) * 100
	# Filter
	identity_threshold=90
	coverage_threshold = 50

	filtered_df = df[(df.iloc[:, 6] <= 200) | (df.iloc[:, 7] <= 200) | (df.iloc[:, 8] <= 200) | (df.iloc[:, 9] <= 200) | 
						(df.iloc[:, 6] >= df.loc[:, "contig_length"]-200) | (df.iloc[:, 7] >= df.loc[:, "contig_length"]-200)]

	filtered_df2 = filtered_df[((filtered_df.iloc[:, 2] >= identity_threshold) & 
								(filtered_df.loc[:, 'total_coverage'] >= coverage_threshold)) | 
							((filtered_df.iloc[:, 6] <= 200) & (filtered_df.iloc[:, 2] >= identity_threshold)) |
							((filtered_df.iloc[:, 7] <= 200) & (filtered_df.iloc[:, 2] >= identity_threshold)) |
							((filtered_df.iloc[:, 6] >= filtered_df.loc[:, "contig_length"]-200) & 
								(filtered_df.iloc[:, 2] >= identity_threshold)) | 
							((filtered_df.iloc[:, 7] >= filtered_df.loc[:, "contig_length"]-200) & 
								(filtered_df.iloc[:, 2] >= identity_threshold))]
	filtered_df2 = filtered_df2.reset_index(drop=True)
	#print("First time")
	#print(filtered_df2)
	# Assuming df4.iloc[:, 2] contains the single string as a list
	rna_list_string = df4.iloc[:, 2].tolist()[0]

	# Split the string into a list of individual items
	rna_list = rna_list_string.split(', ')

	rna_dict = {}

	for item in rna_list:
	    # Split the string to separate the part outside and inside the parentheses
		value, full_id = item.split('(')

        # Remove the closing parenthesis from the full_id and strip any extra whitespace
		full_id = full_id.rstrip(')')

		if " |" in full_id:
			value1, value2 = full_id.split(" |")
		else:
		    value1 = full_id
		    value2 = None  # or any other default value if no separator is present

        # Add to the dictionary
		rna_dict[value1] = item
	
	
	# File path (replace with your file path if needed)
	file_path = db_path+"/genomes/filtered_" + df4.iloc[0, 1] + "_rna_from_genomic.fna"
	#print("File_path : ",file_path)
	# Initialize a list to hold the extracted data
	data = []

	# Variables to hold the current sequence information
	current_id = ''
	current_gene = ''
	current_product = ''
	current_location = ''
	current_sequence = ''

	# Read the file line by line and extract the required information
	with open(file_path, 'r') as file:
		for line in file:
		    if line.startswith('>'):  # Only process header lines starting with '>'
		        # If there is a current sequence, save the previous entry
		        if current_id and current_sequence:
		            # Append the extracted data to the list
		            if current_product in ["5S ribosomal RNA", "16S ribosomal RNA", "23S ribosomal RNA",
		            						"5S Ribosomal RNA", "16S Ribosomal RNA", "23S Ribosomal RNA"]:
		                # Extract component after "_rrna_"
		                gene_id = f"{current_product}({current_id})"
		                data.append([current_id, current_gene, current_product, current_location, current_sequence, gene_id])

		        # Reset sequence for the new record
		        current_sequence = ''
		        
		        # Extract ID
		        id_match = re.search(r'lcl\|(.+?) ', line)
		        current_id = id_match.group(1) if id_match else ''
		        
		        # Extract gene
		        gene_match = re.search(r'\[gene=(.+?)\]', line)
		        current_gene = gene_match.group(1) if gene_match else ''
		        
		        # Extract product
		        product_match = re.search(r'\[product=(.+?)\]', line)
		        current_product = product_match.group(1) if product_match else ''
		        
		        # Extract location
		        location_match = re.search(r'\[location=(.+?)\]', line)
		        current_location = location_match.group(1) if location_match else ''
		        
		    else:
		        # Append the sequence line to the current sequence
		        current_sequence += line.strip()

		# Save the last record if it exists and matches the product filter
		if current_id and current_sequence:
		    if current_product in ["5S ribosomal RNA", "16S ribosomal RNA", "23S ribosomal RNA",
		            						"5S Ribosomal RNA", "16S Ribosomal RNA", "23S Ribosomal RNA"]:
		        # Create gene_id without space between product and the component
		        gene_id = f"{current_product}({current_id})"
		        data.append([current_id, current_gene, current_product, current_location, current_sequence, gene_id])

	# Create a DataFrame with the extracted data
	df = pd.DataFrame(data, columns=['ID', 'Gene', 'Product', 'Location', 'Sequence', 'Gene_ID'])

	# Display the DataFrame
	#print(df)
	
	# Extract the 'Gene_ID' column and convert it to a list
	files = df['Gene_ID'].tolist()

	# Display the list of gene_ids
	#print(files)

	# Initialize an empty dictionary to store the filtered DataFrames
	blastn_data = {}
	#a fisrt column to show how many copies of rRNA genes are in this genome 
	
	# Initialize a dictionary to hold the data
	rRNA_dict = {'Gene': [], 'Number of Copies': [], 'Names of Copies': []}

	# Process the keys to extract information
	rRNA_types = ['5S', '16S', '23S']
	for rRNA_type in rRNA_types:
		# Filter the rRNA genes for the current rRNA type
		rRNA_names = [
				key for key in rna_dict.values()
				if key.startswith(f'{rRNA_type} ribosomal RNA') or key.startswith(f'{rRNA_type} Ribosomal RNA')
			]
		# Add the rRNA type to the 'Gene' column
		rRNA_dict['Gene'].append(f'{rRNA_type} rRNA')
		
		# Add the number of copies to the 'Number of Copies' column
		rRNA_dict['Number of Copies'].append(len(rRNA_names))
		
		# Add the names of the copies to the 'Names of Copies' column
		rRNA_dict['Names of Copies'].append(', '.join(rRNA_names))

	# Create a DataFrame
	rRNA_df = pd.DataFrame(rRNA_dict)

	# Display the DataFrame
	#print("rRNA_df")
	#print(rRNA_df)

	for file in files:
		# The 'file' name itself represents the key we're interested in, so we can use it directly
		rRNA_gene_value = file
		#print(rRNA_gene_value)
		# Filter 'all_results_df' to get only rows where 'rRNA_gene' matches the current 'file' name
		filtered_df2.iloc[:, 1] = filtered_df2.iloc[:, 1].str.replace('lcl\|', '', regex=True)
		filtered_df = filtered_df2[filtered_df2.iloc[:, 16] == rRNA_gene_value]
		
		# Save the filtered DataFrame in the dictionary with the key being the 'file' name
		try:
			blastn_data[filtered_df['rRNA_gene'].iloc[0]] = filtered_df
		except IndexError:
			# Skip adding the entry to the dictionary if there are no hits for the current rRNA gene
			pass

	#print("filtered_df2")
	#print(filtered_df2.iloc[:, 16])
	# Initialize an empty dictionary to store the valid hits for each contig
	valid_hits = {}

	for gene, df in blastn_data.items():
		# Sort the DataFrame by bitscore in descending order
		df = df.sort_values(by=11, ascending=False)
		
		# Iterate over each row in the DataFrame
		for _, row in df.iterrows():
			contig = row[0]
			
			# If the contig is not in the valid_hits dictionary, initialize an empty list for it
			if contig not in valid_hits:
				valid_hits[contig] = []
			
			# Check if the current hit overlaps with any of the previous valid hits for the same contig
			overlaps = False
			for i, valid_hit in enumerate(valid_hits[contig]):
				if (row[6] >= valid_hit[6] and row[6] <= valid_hit[7]) or (row[7] >= valid_hit[6] and row[7] <= valid_hit[7]):
					overlaps = True
					# If the current hit has a better bitscore, replace the overlapping hit with the current hit
					if row[11] > valid_hit[11]:
						valid_hits[contig][i] = row
					break
			
			# If the current hit doesn't overlap with any previous valid hits for the same contig, add it to the valid hits list
			if not overlaps:
				valid_hits[contig].append(row)

	#print("valid_hits")
	#print(valid_hits)
	# Create a new DataFrame to store the valid hits
	data = pd.DataFrame(columns=['contig', 'gene', 'bitscore', 'length', 'sequence', 'start', 'end', 'reference_sequence', 'ref_start', 'ref_end', 'identity', 'alignment_pct', 'e-value', 'gaps', 'mismatches', 'Reference_genome'])

	# Iterate over each contig and its valid hits
	for contig, hits in valid_hits.items():
		for hit in hits:
			# Create a new DataFrame with the extracted information
			new_row = pd.DataFrame({
				'contig': [hit[0]],
				'gene': [hit["rRNA_gene"]],
				'bitscore': [hit[11]],
				'length': [hit[14]],
				'sequence': [hit[12]],
				'start': [hit[6]],
				'end': [hit[7]],
				'reference_sequence': [hit[13]],
				'ref_start': [hit[8]],
				'ref_end': [hit[9]],
				'identity': [hit[2]],
				'alignment_pct': [hit[3] / hit[14] * 100],
				'e-value': [hit[10]],
				'gaps': [hit[5]],
				'mismatches': [hit[4]],
				'Reference_genome': [hit["Reference_Genome"]]
			})
			
			# Concatenate the new row with the existing DataFrame
			data = pd.concat([data, new_row], ignore_index=True)

	#print(data)
	data = data.sort_values(by='gene', ascending=True)
	data = data.dropna(subset=['reference_sequence'])
	final_seqs = ''
	for index, row in data.iterrows():
		final_seqs = str(final_seqs) + str(row['sequence'])
	#print("Print data")
	#print(data)
	
	# Function to correct the strand of the sequence based on the reference sequence 
	
	def reorient_sequences(row):
		sequence = row.iloc[4]

		if row.iloc[8] > row.iloc[9]:  # reference is on reverse strand
		    sequence = str(Seq(sequence).reverse_complement())
		    strand = "negative"
		else:
		    strand = "positive"

		row.iloc[4] = sequence
		row['strand'] = strand
	return row

	data = data.apply(reorient_sequences, axis=1)

	# Function to align sequences and find differences
	# Function to align sequences and find differences
	def find_sequence_differences(row):
		reference_sequence = row['reference_sequence']
		sequence = row['sequence']
		

		if (row.iloc[9] < row.iloc[8] and row.iloc[6] > row.iloc[5]):
					# Reverse complement the sequence using Bio.Seq
			reference_sequence = str(Seq(reference_sequence).reverse_complement())

		if (row.iloc[9] > row.iloc[8] and row.iloc[6] < row.iloc[5]):
			sequence = str(Seq(sequence).reverse_complement())

		aligner = Align.PairwiseAligner()
		aligner.mode = 'global'  # Needleman-Wunsch for global alignment
		aligner.match_score = 1
		aligner.mismatch_score = -1
		aligner.open_gap_score = -2
		aligner.extend_gap_score = -0.5
		seq1= sequence
		seq2= reference_sequence
		alignments = aligner.align(seq1, seq2)
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
		#
		differences = []
		prev_event = None
		prev_position = None
		start_position = None
		continuous_change = ""
		change = ""  # Initialize change variable with an empty string

		for i, (q, s) in enumerate(zip(query_aligned, subject_aligned)):
			if q != s:
				if (q == "-") and (s != "-"):
					event = "deletion"
					change = s  # Update change variable for deletion event
				elif (q != "-") and (s == "-"):
					event = "insertion"
					change = q  # Update change variable for insertion event
				elif (q == "*") and (s != "*"):
					differences.append(f"At position {i + 1}: STOP CODON")
					break
				else:
					event = "substitution"
					change = f"{s} -> {q}"  # Update change variable for substitution event

				if event == prev_event and i == prev_position + 1:
					# Continuous event, update the change
					if event != "substitution":
						continuous_change += change
				else:
					# New event
					if prev_event is not None and prev_event != "substitution":
						# Report the previous continuous event as a range with nucleotides
						end_position = prev_position + 1
						differences[-1] = f"At positions {start_position}-{end_position}: {prev_event} of {continuous_change}"
						continuous_change = ""
					elif prev_event == "substitution":
						differences[-1] = f"At position {prev_position + 1}: {prev_event} of {differences[-1].split(': ')[-1]}"
					
					# Start a new event
					start_position = i + 1
					if event != "substitution":
						continuous_change = change
						differences.append(f"At position {i + 1}: {event}")
					else:
						differences.append(f"At position {i + 1}: {change}")

				prev_event = event
				prev_position = i

		# Report the last continuous event as a range with nucleotides if applicable
		if prev_event is not None and prev_event != "substitution":
			end_position = prev_position + 1
			differences[-1] = f"At positions {start_position}-{end_position}: {prev_event} of {continuous_change}"
		elif prev_event == "substitution":
			differences[-1] = f"At position {prev_position + 1}: {prev_event} of {differences[-1].split(': ')[-1]}"

		value_for_new_column = differences if differences else []
		wrapped_lists = [[lst] for lst in value_for_new_column]
		
		return wrapped_lists if wrapped_lists else 'No changes'
		

	# Apply the function to each row
	data['sequence_changes'] = data.apply(find_sequence_differences, axis=1)
	
	# Function to extract positions and match mutations
	def find_mutations(row, mutations_df, difference_number):
		positions = []
		# Ensure that 'differences' is a list of lists
		if isinstance(row['sequence_changes'], list):
		    for diff in row['sequence_changes']:
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
		        # Extract the part of 'gene' before the parenthesis
		        gene_name = row[1].split('(')[0].strip()
		        # Find all mutations that are within 20 amino acids above or below the position
		        # and have the same Gene as the input row
		        nearby_mutations = mutations_df[
		            (mutations_df['Gene_ID'] == gene_name) &     # Check for matching Gene
		            mutations_df['Codon_pos'].notna() &         # Exclude NaN values
		            (mutations_df['Codon_pos'] != 0) &  
		            mutations_df['Codon_pos'].between(pos - difference_number, pos + difference_number)
		        ]['Mutations_Resistance'].tolist()
		        #print(nearby_mutations)
		        matched_mutations.extend(nearby_mutations)
		        #print(nearby_mutations)
		#return (matched_mutations)
		# Return the matched mutations as a semicolon-separated string, or indicate no matches found
		return "; ".join(matched_mutations) if matched_mutations else "No known mutations found"

	# Function to filter species rows based on exact match or partial match
	def filter_species(row, species_short):
		# Normalize the species name by replacing spaces with underscores
		species_value = str(row['Species']).replace(' ', '_')  # Replace spaces with underscores
		# Check if the species is either exactly or partially matching the normalized Species column
		return species_short in species_value or species_value in species_short

	mutations_details=pd.read_csv(db_path+'full_list_mutations.csv',
    sep=',',                    # Specify the delimiter if it's different from a comma
)

	# Step 1: Convert the 'Values' column to numeric, coercing errors to NaN for non-numeric values
	mutations_details['Codon_pos'] = pd.to_numeric(mutations_details['Codon_pos'], errors='coerce')  # Correct usage on a column, not the entire DataFrame

	# Step 2: Fill NaN values with a default integer (e.g., 0) before converting to int
	mutations_details['Codon_pos'] = mutations_details['Codon_pos'].fillna(0).astype(int)

	species_short = '_'.join(species.split(' ')[:2])


	# Remove rows where 'Species' contains 'Mycobacterium'
	mutations_details = mutations_details[~mutations_details['Species'].str.contains('Mycobacterium', na=False, case=False)]


	# Normalize the Species column in the DataFrame to replace spaces with underscores
	mutations_details['Normalized_Species'] = mutations_details['Species'].str.replace(' ', '_')

	# Create a mask for rows that match the species criteria
	match_mask = mutations_details['Normalized_Species'].notna() & mutations_details.apply(
		lambda row: filter_species(row, species_short), axis=1
	)

	mutations_details_Species = mutations_details[match_mask]

	# Create the DataFrame with non-matching rows
	mutations_details_OtherSpecies = mutations_details[~match_mask]


		# Step 1: Extract the first two words of the species name and format with an underscore
	formatted_species = "_".join(species.split()[:2])

	# Step 2: Dynamically create the column name based on the formatted species
	column_name = f"{formatted_species} Known Mutations                                                                                                                                                                                                                                                                               ."

	# Apply the function to create the 'Species_Mutations' column
	data[column_name] = data.apply(
		lambda row: find_mutations(row, mutations_details_Species,difference_number), axis=1
	)

	# Apply the function to create the 'Species_Mutations' column
	data['Other_Species_Known_Mutations']=data.apply(
		lambda row: find_mutations(row, mutations_details_OtherSpecies,difference_number), axis=1
	)

	# Step 1: Get all rRNA gene names from rRNA_df
	#all_rRNA_copies = []
	#for names in rRNA_df['Names of Copies']:
	#	all_rRNA_copies.extend(names.split(', '))

	# Step 2: Check which rRNA copies are in gene_data_df
	#existing_copies = data['gene'].tolist()

	# Step 3: Find the missing copies
	#missing_copies = [copy for copy in all_rRNA_copies if copy not in existing_copies]

	# Step 4: Add missing copies to gene_data_df, marking them as "MISSING COPY"
	#for missing_copy in missing_copies:
	#	missing_row = {
	#	    'contig': 'MISSING COPY',
	#	    'gene': missing_copy,
	#	}
	#	data = pd.concat([data, pd.DataFrame([missing_row])], ignore_index=True)

	# Step 5: Fill other columns with "MISSING COPY" for these new rows
	#for col in data.columns:
	#	data[col].fillna('MISSING COPY', inplace=True)
	data = data.sort_values(by='gene', ascending=True)
	
	data = data.loc[:,['contig', 'gene', "sequence_changes",column_name, 
					"Other_Species_Known_Mutations", "sequence", "reference_sequence",
					"start","end", "identity", "alignment_pct", "length","e-value", "bitscore",
					"gaps", "mismatches", "Reference_genome"]]
 					
	data.rename(columns={"contig":"Contig ID", 'gene': 'rRNA gene ID', 
				"sequence_changes": "Differences from reference gene                                                                                                            .",
				"start" : "Query Start", "end" : "Query Stop",
    			"identity" : "Identity", 
    			"Other_Species_Known_Mutations" : "Other Species Known Mutations                                                                                                            .",
				"alignment_pct" :  "Coverage",
				"e-value" : "E-value",
				"bitscore" : "BitScore",
				"sequence" : "rRNA Copy Sequence                                                                                                                                                                                                                                                                                                                                                                                                         .",
    			"reference_sequence" : "rRNA Copy Reference Sequence                                                                                                                                                                                                                                                                                                                                                                                                         .",
				"length" : "Reference rRNA gene length",
				"gaps" : "Gaps" , "mismatches" : "Mismatches"}, inplace=True) 
    			
    			# Specify the column you want to keep in its original data type
	column_to_exclude = ["Identity","Coverage",
    				 "E-value","BitScore"]


    # Convert all columns to strings without decimals, except the specified columns
	data = data.apply(lambda x: x.astype(int).astype(str) if x.name not in column_to_exclude and x.dtype != 'object' else x)
	data["Coverage"] = pd.to_numeric(data["Coverage"], errors='coerce')
	data.loc[data["Coverage"] > 100, "Coverage"] = 100

	# Display the updated gene_data_df
	return(data)
