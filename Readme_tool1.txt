 Here are the scripts for vas to run the 1st component of the tool in python 
 
 The folder databases has all the needed files. 
 
 In python run : 

query_file = "S3323_contigs.fasta"  # Replace with the actual fasta file from the submit
db_path = "databases/all/all_amr"
num_threads = 4
protein_annotation_file="databases/genes_annotation_databases.csv"
protein_start_filter = 20 # Threshold for protein start position filter.
identity_threshold= 70 # Threshold for sequence identity.
coverage_threshold= 70 # Threshold for coverage percentage.



import find_amr_genes_module

blast_results_df = find_amr_genes_module.run_blastx_to_dataframe(query_file, db_path, num_threads)

final_results=find_amr_genes_module.process_blast_results(protein_annotation_file, blast_results_df, query_file, protein_start_filter, identity_threshold, coverage_threshold)


