For finding the species for the tool are 
cat refseq/genomes_assemblies_CoreGenes_final.txt | awk -F"\t" '{print $1}' | sort > refseq/genomes_assemblies_CoreGenes_final_species.txt 
cat db/genomes_assemblies_rRNA_final.txt | awk -F"\t" '{print $1}'| sort > db/genomes_assemblies_rRNA_species.txt
 comm -12 refseq/genomes_assemblies_CoreGenes_final_species.txt db/genomes_assemblies_rRNA_species.txt > available_species.txt

In python 

#Tool 1 

import sys
import importlib.util

# Define the module path
module_path = '/home/argis/Desktop/pasteur/amr_server/AmrProfiler_new/find_amr_genes_module.py'

# Load the module by specifying the module name as a string
spec = importlib.util.spec_from_file_location("find_amr_genes_module", module_path)
find_amr_genes_module = importlib.util.module_from_spec(spec)
sys.modules["find_amr_genes_module"] = find_amr_genes_module
spec.loader.exec_module(find_amr_genes_module)


query_file = "21.fasta"
Database_of_Genes = "Resfinder+ReferenceGeneCatalog" #auto einai pou tha dinei o xristis 

if Database_of_Genes == "Resfinder+ReferenceGeneCatalog":
    db_path = "/home/argis/Desktop/pasteur/amr_server/AmrProfiler_new/databases/all/amrFinder_ResFinder"
elif Database_of_Genes == "Resfinder+ReferenceGeneCatalog+Card":  # Replace with the correct alternative condition
    db_path = "/home/argis/Desktop/pasteur/amr_server/AmrProfiler_new/databases/all/all_amr"
#db_path = "/home/argis/Desktop/pasteur/amr_server/AmrProfiler_new/databases/all/all_amr"
num_threads = 4
protein_annotation_file="/home/argis/Desktop/pasteur/amr_server/AmrProfiler_new/databases/genes_annotation_databases.csv"
protein_start_filter = 50 # Threshold for protein start position filter.
identity_threshold= 70 # Threshold for sequence identity.
coverage_threshold= 70 # Threshold for coverage percentage.


blast_results_df = find_amr_genes_module.run_blastx_to_dataframe(query_file, db_path, num_threads)

final_results=find_amr_genes_module.process_blast_results(protein_annotation_file, blast_results_df, query_file, protein_start_filter, identity_threshold, coverage_threshold)

blast_results_df.to_csv("blast_results.csv", index=False, float_format="%.6f")
final_results.to_csv("final_results_tool1.csv", index=False, float_format="%.6f")



#Tool 2 

# Define the module path
module_path = '/home/argis/Desktop/pasteur/amr_server/AmrProfiler/core_genes_finder.py'

# Load the module by specifying the module name as a string
spec = importlib.util.spec_from_file_location("core_genes_finder", module_path)
core_genes_finder = importlib.util.module_from_spec(spec)
sys.modules["core_genes_finder"] = core_genes_finder
spec.loader.exec_module(core_genes_finder)

query_file = "21.fasta"
species = "Pseudomonas aeruginosa"
db_path_folder = "/home/argis/Desktop/pasteur/amr_server/AmrProfiler_new/refseq/"
threads = 4
difference_number = 2  # Example threshold.
difference_from_start = 20
identity_threshold = 90
coverage_threshold = 90


df4=core_genes_finder.get_species_row(species, db_path_folder)  #(the second table that we will show to the user)
df4.to_csv("core_genes_of_this_species.csv", index=False, float_format="%.6f")
blast_results_core=core_genes_finder.run_blastx_for_references(query_file, df4, db_path_folder)
blast_results_core.to_csv("blast_results_core.csv", index=False, float_format="%.6f")

result_df = core_genes_finder.process_blast_results(species, blast_results_core, query_file, df4, difference_number,db_path_folder, difference_from_start,identity_threshold,coverage_threshold)
result_df.to_csv("mutations_results.csv", index=False, float_format="%.6f")




#Tool 3 
module_path = '/home/argis/Desktop/pasteur/amr_server/AmrProfiler/rRNA_genes_finder_NEW.py'

# Load the module by specifying the module name as a string
spec = importlib.util.spec_from_file_location("rRNA_genes_finder_NEW", module_path)
rRNA_genes_finder_NEW = importlib.util.module_from_spec(spec)
sys.modules["rRNA_genes_finder_NEW"] = rRNA_genes_finder_NEW
spec.loader.exec_module(rRNA_genes_finder_NEW)

query_file = "21.fasta"
species = "Pseudomonas aeruginosa"
db_path_folder = "/home/argis/Desktop/pasteur/amr_server/AmrProfiler_new/db/"
difference_number_rRNA = 5 


df5=rRNA_genes_finder_NEW.get_species_row(species,db_path_folder)  #(the second table that we will show to the user)
df5.to_csv("rRNA_genes_of_this_species.csv", index=False, float_format="%.6f")

blast_results_rRNA=rRNA_genes_finder_NEW.run_blastn_for_references(query_file, df5,db_path_folder)
rRNA_dif = rRNA_genes_finder_NEW.find_rRNA_differences(species,df5, blast_results_rRNA, query_file, db_path_folder, difference_number_rRNA,5)
rRNA_dif.to_csv("mutations_rRNA_results.csv", index=False, float_format="%.6f")

#To run automatically with default parameters 
python /home/argis/Desktop/pasteur/amr_server/AmrProfiler_new/amrprofiler.py "Lar4933_contigs.fasta" "Staphylococcus aureus" /home/argis/Desktop/pasteur/amr_server/AmrProfiler_new/
