#Download refseq on 04/07/2024 
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt

awk -F "\t" '($5 == "representative genome" || $5 == "reference genome") {print $8"\t"$20}' assembly_summary.txt > ftpdirpaths.txt

#the suffix needed is _rna_from_genomic.fna.gz 
awk 'BEGIN{FS=OFS="/";filesuffix="rna_from_genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' ftpdirpaths.txt > genomes_to_download.txt

cat genomes_to_download.txt | awk -F"\t" '{print $2}' > for_download.txt 

mkdir genomes 

bash download_rRNA.sh 

bash create_blast_databases.sh

We have to copy the full_list_mutations.csv from the and ran the script process_ribosomal_RNAs_NEW.ipynb
/home/argis/Desktop/pasteur/amr_server/AmrProfiler/refseq/all_detailed to the folder db
and then run 

df4=get_species_row(species)  #(the second table that we will show to the user)

blast_results=run_blastn_for_references(query_file, df4)

def find_rRNA_differences(df4, blast_results, query_file, reference_df, num_threads=5):
