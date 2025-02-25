#Copy the file genomes_assemblies_CoreGenes.txt from the folder refseq

#Run in python 

species = "Enterococcus faecalis"
query_file = "GCF_029024925.1_ASM2902492v1_genomic.fna"
difference_number = 2  # Example threshold.
difference_from_start = 20
identity_threshold = 70
coverage_threshold = 70

import core_genes_finder

df4=core_genes_finder.get_species_row(species)  #(the second table that we will show to the user)

blast_results=core_genes_finder.run_blastx_for_references(query_file, df4)
result_df = core_genes_finder.process_blast_results(species, blast_results, query_file, df4, difference_number,difference_from_start,identity_threshold,coverage_threshold)
#(the third table that we will show to the user)


