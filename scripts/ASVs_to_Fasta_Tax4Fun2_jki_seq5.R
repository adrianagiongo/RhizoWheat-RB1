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

# Function to read the CSV file containing ASV names to be selected
read_asv_names <- function(csv_file) {
  # Read the CSV file
  df <- read.csv(csv_file)
  
  # Ensure the dataframe has the correct column
  if (!"ASV" %in% colnames(df)) {
    stop("The CSV file must contain an 'ASV' column.")
  }
  
  # Extract ASV names
  asv_names <- df$ASV
  
  return(asv_names)
}

# Function to save selected sequences as a fasta file
save_selected_asvs_to_fasta <- function(sequences, selected_names, output_file) {
  # Filter sequences based on the provided list of names
  selected_sequences <- sequences[selected_names]
  
  # Remove NA values
  selected_sequences <- selected_sequences[!is.na(selected_sequences)]
  
  # Print the selected sequences for debugging
  print("Selected sequences:")
  print(selected_sequences)
  
  # Create DNAStringSet object
  dna_string_set <- DNAStringSet(selected_sequences)
  names(dna_string_set) <- names(selected_sequences)
  
  # Write to fasta file
  writeXStringSet(dna_string_set, filepath = output_file)
}


#### T1_D30
# Specify the path to the sequence file
csv_file_sequences <- "~/Documents/R_analysis/jki_seq5/data_jki_seq5/Tax4Fun2/jki_seq5_seqs_Tax4Fun.csv"  # Path to your CSV file with sequences

# Specify the path to the file containing ASV names to be selected
csv_file_asv_names <- "~/Documents/R_analysis/jki_seq5/data_jki_seq5/Tax4Fun2/T1_D30_filt_deseq2.csv"  # Path to your CSV file with ASV names

# Read the sequences from the CSV file
sequences <- read_sequences_from_csv(csv_file_sequences)

# Print the loaded sequences for debugging
print("Loaded sequences:")
print(sequences)

# Read the ASV names to select
asv_names_to_select <- read_asv_names(csv_file_asv_names)

# Print the ASV names to select for debugging
print("ASV names to select:")
print(asv_names_to_select)

# Filter out ASV names that do not exist in the sequences
existing_asv_names <- asv_names_to_select[asv_names_to_select %in% names(sequences)]
missing_asv_names <- setdiff(asv_names_to_select, names(sequences))

# Print the existing and missing ASV names for debugging
print("Existing ASV names:")
print(existing_asv_names)
print("Missing ASV names:")
print(missing_asv_names)

# Print warning for any missing ASV names
if (length(missing_asv_names) > 0) {
  warning("The following ASV names were not found in the sequences: ", paste(missing_asv_names, collapse = ", "))
}

# Specify the output fasta file name
output_fasta_file <- "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/T1_D30_filt_deseq2.fasta"

# Save the selected sequences to a fasta file
save_selected_asvs_to_fasta(sequences, existing_asv_names, output_fasta_file)

###############
#### T1_D60
# Specify the path to the sequence file
csv_file_sequences <- "~/Documents/R_analysis/jki_seq5/data_jki_seq5/Tax4Fun2/jki_seq5_seqs_Tax4Fun.csv"  # Path to your CSV file with sequences

# Specify the path to the file containing ASV names to be selected
csv_file_asv_names <- "~/Documents/R_analysis/jki_seq5/data_jki_seq5/Tax4Fun2/T1_D60_filt_deseq2.csv"  # Path to your CSV file with ASV names

# Read the sequences from the CSV file
sequences <- read_sequences_from_csv(csv_file_sequences)

# Print the loaded sequences for debugging
print("Loaded sequences:")
print(sequences)

# Read the ASV names to select
asv_names_to_select <- read_asv_names(csv_file_asv_names)

# Print the ASV names to select for debugging
print("ASV names to select:")
print(asv_names_to_select)

# Filter out ASV names that do not exist in the sequences
existing_asv_names <- asv_names_to_select[asv_names_to_select %in% names(sequences)]
missing_asv_names <- setdiff(asv_names_to_select, names(sequences))

# Print the existing and missing ASV names for debugging
print("Existing ASV names:")
print(existing_asv_names)
print("Missing ASV names:")
print(missing_asv_names)

# Print warning for any missing ASV names
if (length(missing_asv_names) > 0) {
  warning("The following ASV names were not found in the sequences: ", paste(missing_asv_names, collapse = ", "))
}

# Specify the output fasta file name
output_fasta_file <- "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/T1_D60_filt_deseq2.fasta"

# Save the selected sequences to a fasta file
save_selected_asvs_to_fasta(sequences, existing_asv_names, output_fasta_file)

##############

#### T2_D30
# Specify the path to the sequence file
csv_file_sequences <- "~/Documents/R_analysis/jki_seq5/data_jki_seq5/Tax4Fun2/jki_seq5_seqs_Tax4Fun.csv"  # Path to your CSV file with sequences

# Specify the path to the file containing ASV names to be selected
csv_file_asv_names <- "~/Documents/R_analysis/jki_seq5/data_jki_seq5/Tax4Fun2/T2_D30_filt_deseq2.csv"  # Path to your CSV file with ASV names

# Read the sequences from the CSV file
sequences <- read_sequences_from_csv(csv_file_sequences)

# Print the loaded sequences for debugging
print("Loaded sequences:")
print(sequences)

# Read the ASV names to select
asv_names_to_select <- read_asv_names(csv_file_asv_names)

# Print the ASV names to select for debugging
print("ASV names to select:")
print(asv_names_to_select)

# Filter out ASV names that do not exist in the sequences
existing_asv_names <- asv_names_to_select[asv_names_to_select %in% names(sequences)]
missing_asv_names <- setdiff(asv_names_to_select, names(sequences))

# Print the existing and missing ASV names for debugging
print("Existing ASV names:")
print(existing_asv_names)
print("Missing ASV names:")
print(missing_asv_names)

# Print warning for any missing ASV names
if (length(missing_asv_names) > 0) {
  warning("The following ASV names were not found in the sequences: ", paste(missing_asv_names, collapse = ", "))
}

# Specify the output fasta file name
output_fasta_file <- "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/T2_D30_filt_deseq2.fasta"

# Save the selected sequences to a fasta file
save_selected_asvs_to_fasta(sequences, existing_asv_names, output_fasta_file)

###############
#### T2_D60
# Specify the path to the sequence file
csv_file_sequences <- "~/Documents/R_analysis/jki_seq5/data_jki_seq5/Tax4Fun2/jki_seq5_seqs_Tax4Fun.csv"  # Path to your CSV file with sequences

# Specify the path to the file containing ASV names to be selected
csv_file_asv_names <- "~/Documents/R_analysis/jki_seq5/data_jki_seq5/Tax4Fun2/T2_D60_filt_deseq2.csv"  # Path to your CSV file with ASV names

# Read the sequences from the CSV file
sequences <- read_sequences_from_csv(csv_file_sequences)

# Print the loaded sequences for debugging
print("Loaded sequences:")
print(sequences)

# Read the ASV names to select
asv_names_to_select <- read_asv_names(csv_file_asv_names)

# Print the ASV names to select for debugging
print("ASV names to select:")
print(asv_names_to_select)

# Filter out ASV names that do not exist in the sequences
existing_asv_names <- asv_names_to_select[asv_names_to_select %in% names(sequences)]
missing_asv_names <- setdiff(asv_names_to_select, names(sequences))

# Print the existing and missing ASV names for debugging
print("Existing ASV names:")
print(existing_asv_names)
print("Missing ASV names:")
print(missing_asv_names)

# Print warning for any missing ASV names
if (length(missing_asv_names) > 0) {
  warning("The following ASV names were not found in the sequences: ", paste(missing_asv_names, collapse = ", "))
}

# Specify the output fasta file name
output_fasta_file <- "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/T2_D60_filt_deseq2.fasta"

# Save the selected sequences to a fasta file
save_selected_asvs_to_fasta(sequences, existing_asv_names, output_fasta_file)

