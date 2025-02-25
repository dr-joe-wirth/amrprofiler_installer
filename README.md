# AmrProfiler

**AmrProfiler** is a bioinformatics tool designed to identify antimicrobial resistance (AMR) genes, detect point mutations in core genes, and analyze mutations in rRNA genes across 18,000 bacterial species.

## Installation & Setup

To run AmrProfiler locally, you must first download the required databases and files from [Zenodo](#) (replace with actual link). Extract the downloaded folders inside the `amrprofiler` directory.

## Running AmrProfiler

### Default Execution

To run AmrProfiler with default parameters, execute:

```sh
python /home/argis/Desktop/pasteur/amr_server/AmrProfiler_new/amrprofiler.py "Lar4933_contigs.fasta" "Staphylococcus aureus"
```

### Tool 1: AMR Gene Identification

This tool identifies antimicrobial resistance genes from a given query sequence.

```python
import sys
import importlib.util

# Load module
module_path = '/home/argis/Desktop/pasteur/amr_server/AmrProfiler_new/find_amr_genes_module.py'
spec = importlib.util.spec_from_file_location("find_amr_genes_module", module_path)
find_amr_genes_module = importlib.util.module_from_spec(spec)
sys.modules["find_amr_genes_module"] = find_amr_genes_module
spec.loader.exec_module(find_amr_genes_module)

# Parameters
query_file = "21.fasta"
Database_of_Genes = "Resfinder+ReferenceGeneCatalog"
db_path = "/home/argis/Desktop/pasteur/amr_server/AmrProfiler_new/databases/all/amrFinder_ResFinder" if Database_of_Genes == "Resfinder+ReferenceGeneCatalog" else "/home/argis/Desktop/pasteur/amr_server/AmrProfiler_new/databases/all/all_amr"
num_threads = 4
protein_annotation_file = "/home/argis/Desktop/pasteur/amr_server/AmrProfiler_new/databases/genes_annotation_databases.csv"
protein_start_filter = 50
identity_threshold = 70
coverage_threshold = 70

# Run analysis
blast_results_df = find_amr_genes_module.run_blastx_to_dataframe(query_file, db_path, num_threads)
final_results = find_amr_genes_module.process_blast_results(protein_annotation_file, blast_results_df, query_file, protein_start_filter, identity_threshold, coverage_threshold)

# Save results
blast_results_df.to_csv("blast_results.csv", index=False, float_format="%.6f")
final_results.to_csv("final_results_tool1.csv", index=False, float_format="%.6f")
```

### Tool 2: Core Gene Mutation Analysis

This tool detects point mutations in core genes for a given bacterial species.

```python
# Load module
module_path = '/home/argis/Desktop/pasteur/amr_server/AmrProfiler/core_genes_finder.py'
spec = importlib.util.spec_from_file_location("core_genes_finder", module_path)
core_genes_finder = importlib.util.module_from_spec(spec)
sys.modules["core_genes_finder"] = core_genes_finder
spec.loader.exec_module(core_genes_finder)

# Parameters
query_file = "21.fasta"
species = "Pseudomonas aeruginosa"
db_path_folder = "/home/argis/Desktop/pasteur/amr_server/AmrProfiler_new/refseq/"
threads = 4
difference_number = 2
difference_from_start = 20
identity_threshold = 90
coverage_threshold = 90

# Run analysis
df4 = core_genes_finder.get_species_row(species, db_path_folder)
df4.to_csv("core_genes_of_this_species.csv", index=False, float_format="%.6f")
blast_results_core = core_genes_finder.run_blastx_for_references(query_file, df4, db_path_folder)
blast_results_core.to_csv("blast_results_core.csv", index=False, float_format="%.6f")

result_df = core_genes_finder.process_blast_results(species, blast_results_core, query_file, df4, difference_number, db_path_folder, difference_from_start, identity_threshold, coverage_threshold)
result_df.to_csv("mutations_results.csv", index=False, float_format="%.6f")
```

### Tool 3: rRNA Mutation Analysis

This tool identifies mutations in rRNA genes for a given bacterial species.

```python
# Load module
module_path = '/home/argis/Desktop/pasteur/amr_server/AmrProfiler/rRNA_genes_finder_NEW.py'
spec = importlib.util.spec_from_file_location("rRNA_genes_finder_NEW", module_path)
rRNA_genes_finder_NEW = importlib.util.module_from_spec(spec)
sys.modules["rRNA_genes_finder_NEW"] = rRNA_genes_finder_NEW
spec.loader.exec_module(rRNA_genes_finder_NEW)

# Parameters
query_file = "21.fasta"
species = "Pseudomonas aeruginosa"
db_path_folder = "/home/argis/Desktop/pasteur/amr_server/AmrProfiler_new/db/"
difference_number_rRNA = 5

# Run analysis
df5 = rRNA_genes_finder_NEW.get_species_row(species, db_path_folder)
df5.to_csv("rRNA_genes_of_this_species.csv", index=False, float_format="%.6f")
blast_results_rRNA = rRNA_genes_finder_NEW.run_blastn_for_references(query_file, df5, db_path_folder)
rRNA_dif = rRNA_genes_finder_NEW.find_rRNA_differences(species, df5, blast_results_rRNA, query_file, df5, db_path_folder, difference_number_rRNA, 5)
rRNA_dif.to_csv("mutations_rRNA_results.csv", index=False, float_format="%.6f")
```

## Output Files

- **blast\_results.csv**: Raw BLAST results from AMR gene search
- **final\_results\_tool1.csv**: Processed AMR gene analysis results
- **core\_genes\_of\_this\_species.csv**: Core genes reference for the species
- **blast\_results\_core.csv**: BLAST results for core genes
- **mutations\_results.csv**: Mutations detected in core genes
- **rRNA\_genes\_of\_this\_species.csv**: Reference rRNA genes for the species
- **mutations\_rRNA\_results.csv**: Mutations detected in rRNA genes

## Dependencies

- Python 3.x
- BLAST+
- Pandas
- Biopython

## License

This project is licensed under the MIT License.

---

This README provides a structured and clearer explanation of how to use AmrProfiler. Let me know if you need any refinements!

