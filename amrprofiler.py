#Tool 1 
import sys
import importlib.util
import warnings
import pandas as pd
warnings.filterwarnings("ignore")

query_file = sys.argv[1]
species = sys.argv[2]
amrprofiler_folder = sys.argv[3]

print("Tool 1")

# Add the module's directory to sys.path
sys.path.append(amrprofiler_folder)

# Import the module
import find_amr_genes_module

db_path = amrprofiler_folder+"/databases/all/nr_all_amr"
num_threads = 4
protein_annotation_file=amrprofiler_folder+"/databases/genes_annotation_databases_ASreviewed.csv"
protein_start_filter = 20 # Threshold for protein start position filter.
identity_threshold= 70 # Threshold for sequence identity.
coverage_threshold= 70 # Threshold for coverage percentage.

print("StartBlastx")
try:
    blast_results_df = find_amr_genes_module.run_blastx_to_dataframe(query_file, db_path, num_threads)
    #blast_results_df.to_csv("blast_results_AMR.csv", index=False, float_format="%.6f")
    final_results = find_amr_genes_module.process_blast_results(protein_annotation_file, blast_results_df, query_file, 
                                                              protein_start_filter, identity_threshold, coverage_threshold)
    final_results.to_csv("final_results_tool1.csv", index=False, float_format="%.6f")
except Exception as e:
    print(f"Error in Tool 1: {str(e)}")
    final_results = pd.DataFrame({"Error": ["No acquired AMR genes Found!"]})
    final_results.to_csv("AMR_genes.csv", index=False, float_format="%.6f")



#Tool 2 
print("Tool 2")

import core_genes_finder

db_path_folder = amrprofiler_folder+"/refseq/"
threads = 4
difference_number = 1  # Example threshold.
difference_from_start = 50
identity_threshold = 90
coverage_threshold = 50


try:
    df4 = core_genes_finder.get_species_row(species, db_path_folder)  #(the second table that we will show to the user)
    if df4.empty:
        raise ValueError("Empty species data")
    df4.to_csv("core_genes_of_this_species.csv", index=False, float_format="%.6f")
    blast_results_core = core_genes_finder.run_blastx_for_references(query_file, df4, db_path_folder)
    blast_results_core.to_csv("blast_results_core.csv", index=False, float_format="%.6f")
except Exception as e:
    print("Most probably the species that you provide is not included in the species list (CHECK FOR MISPELLING)")
    df4 = pd.DataFrame({"Error": ["Species not found or invalid"]})
    df4.to_csv("core_genes_of_this_species.csv", index=False, float_format="%.6f")
    blast_results_core = pd.DataFrame()  # Empty DataFrame for downstream handling

try:
    result_df = core_genes_finder.process_blast_results(species, blast_results_core, query_file, df4, 
                                                      difference_number, db_path_folder, difference_from_start,
                                                      identity_threshold, coverage_threshold)
    if result_df.empty:
        raise ValueError("No mutations on core genes of the selected species")
    result_df.to_csv("Core_mutations_results.csv", index=False, float_format="%.6f")
except Exception as e:
    result_df = pd.DataFrame({"Error": ["No matches Found, Probably Wrong Species Selected!"]})
    result_df.to_csv("Core_mutations_results.csv", index=False, float_format="%.6f")




#Tool 3
print("Tool 3")
module_path = '/home/argis/Desktop/pasteur/amr_server/AmrProfiler_new/rRNA_genes_finder_NEW.py'

import rRNA_genes_finder_NEW

db_path_folder = amrprofiler_folder+"/db/"
difference_number_rRNA = 5 


try:
    df5 = rRNA_genes_finder_NEW.get_species_row(species, db_path_folder)  #(the second table that we will show to the user)
    if df5.empty:
        raise ValueError("Empty species data")
    df5.to_csv("rRNA_genes_of_this_species.csv", index=False, float_format="%.6f")
    blast_results_rRNA = rRNA_genes_finder_NEW.run_blastn_for_references(query_file, df5, db_path_folder)
    blast_results_rRNA.to_csv("blast_results_rRNA.csv", index=False, float_format="%.6f")
except Exception as e:
    print("Most probably the species that you provide is not included in the species list (CHECK FOR MISPELLING)")
    df5 = pd.DataFrame({"Error": ["Species not found or invalid"]})
    df5.to_csv("rRNA_genes_of_this_species.csv", index=False, float_format="%.6f")
    blast_results_rRNA = pd.DataFrame()  # Empty DataFrame for downstream handling

try:
    rRNA_dif = rRNA_genes_finder_NEW.find_rRNA_differences(species, df5, blast_results_rRNA, query_file, 
                                                          db_path_folder, difference_number_rRNA, 5)
    if rRNA_dif.empty:
        raise ValueError("No rRNA genes of the selected species found")
    rRNA_dif.to_csv("rRNA_mutations_results.csv", index=False, float_format="%.6f")
except Exception as e:
    print(f"Error in Tool 3: {str(e)}")
    rRNA_dif = pd.DataFrame({"Error": ["No matches Found, Probably Wrong Species Selected!"]})
    rRNA_dif.to_csv("rRNA_mutations_results.csv", index=False, float_format="%.6f")



