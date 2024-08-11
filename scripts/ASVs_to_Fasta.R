# Load necessary libraries
library(Biostrings)
library(dplyr)

# Function to read the CSV file and load sequences
    read_sequences_from_csv <- function(csv_file) {
      # Read the CSV file
      df <- read.csv(csv_file)
      
      # Ensure the dataframe has the correct columns
      if (!all(c("ASV", "Sequence") %in% colnames(df))) {
        stop("The CSV file must contain 'ASV' and 'Sequence' columns.")
      }
      
      # Create a named vector of sequences
      sequences <- setNames(df$Sequence, df$ASV)
      
      return(sequences)
    }
    
# Function to save selected sequences as a fasta file
save_selected_asvs_to_fasta <- function(sequences, selected_names, output_file) {
      # Filter sequences based on the provided list of names
      selected_sequences <- sequences[selected_names]
      
      # Create DNAStringSet object
      dna_string_set <- DNAStringSet(selected_sequences)
      names(dna_string_set) <- selected_names
      
      # Write to fasta file
      writeXStringSet(dna_string_set, filepath=output_file)
    }
    
######Specify the path to sequence file
csv_file <- "~/Documents/R_analysis/jki_seq5/data_jki_seq5/jki_seq5_seqs.csv"   # Path to your CSV file
    
#Read the sequences from the CSV file
sequences <- read_sequences_from_csv(csv_file)
    
#List of ASV names to select
asv_names_to_select <- c("ASV62", "ASV112", "ASV159",
                          "ASV460", "ASV849", "ASV936",
                          "ASV1106", "ASV1161", "ASV2467",
                          "ASV3591", "ASV6585", "ASV7048",
                          "ASV9978", "ASV10621", "ASV11951",
                          "ASV12221", "ASV12236")  # Update this list with your ASV names
    
#Specify the output fasta file name
output_fasta_file <- "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/selected_asvs.fasta"
    
#Save the selected sequences to a fasta file
save_selected_asvs_to_fasta(sequences, asv_names_to_select, output_fasta_file)
    
    