##Source
#https://github.com/songweizhi/Tax4Fun2_short_tutorial?tab=readme-ov-file

##Making functional predictions using the default database. --> using results from DESeq2

library("Tax4Fun2")

## T1_D30_deseq2
# modify the following 8 lines as needed
query_otu_seq   = '~/Documents/R_analysis/jki_seq5/data_jki_seq5/Tax4Fun2/T1_D30_filt_deseq2.fasta'              # input OTU sequence file
query_otu_table = '~/Documents/R_analysis/jki_seq5/data_jki_seq5/Tax4Fun2/T1_D30_filt_deseq2.txt'          # input OTU table
pwd_op_folder   = '~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tax4Fun2_jki_seq5/T1_D30_deseq2'     # output directory
pwd_ref_data    = '~/Documents/R_analysis/jki_seq5/data_jki_seq5/Tax4Fun2/Tax4Fun2_ReferenceData_v2/'  # path to the default database (i.e., Tax4Fun2_ReferenceData_v2), need to be decompressed before use
norm_by_cn      = TRUE                          # normalize_by_copy_number (TRUE or FALSE)
norm_path       = TRUE                          # normalize_pathways (TRUE or FALSE)
iden            = 0.97                          # min_identity_to_reference, please modify as needed
num_of_threads  = 1                             # number of CPU cores to use, please modify as needed

# predict functions
runRefBlast(path_to_otus = query_otu_seq, path_to_reference_data = pwd_ref_data, path_to_temp_folder = pwd_op_folder, database_mode = "Ref99NR", use_force = T, num_threads = num_of_threads)

makeFunctionalPrediction(path_to_otu_table = query_otu_table, path_to_reference_data = pwd_ref_data, path_to_temp_folder = pwd_op_folder, database_mode = "Ref99NR", normalize_by_copy_number = norm_by_cn, min_identity_to_reference = iden, normalize_pathways = norm_path)

## T1_D60_deseq2
# modify the following 8 lines as needed
query_otu_seq   = '~/Documents/R_analysis/jki_seq5/data_jki_seq5/Tax4Fun2/T1_D60_filt_deseq2.fasta'              # input OTU sequence file
query_otu_table = '~/Documents/R_analysis/jki_seq5/data_jki_seq5/Tax4Fun2/T1_D60_filt_deseq2.txt'          # input OTU table
pwd_op_folder   = '~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tax4Fun2_jki_seq5/T1_D60_deseq2'     # output directory
pwd_ref_data    = '~/Documents/R_analysis/jki_seq5/data_jki_seq5/Tax4Fun2/Tax4Fun2_ReferenceData_v2/'  # path to the default database (i.e., Tax4Fun2_ReferenceData_v2), need to be decompressed before use
norm_by_cn      = TRUE                          # normalize_by_copy_number (TRUE or FALSE)
norm_path       = TRUE                          # normalize_pathways (TRUE or FALSE)
iden            = 0.97                          # min_identity_to_reference, please modify as needed
num_of_threads  = 1                             # number of CPU cores to use, please modify as needed

# predict functions
runRefBlast(path_to_otus = query_otu_seq, path_to_reference_data = pwd_ref_data, path_to_temp_folder = pwd_op_folder, database_mode = "Ref99NR", use_force = T, num_threads = num_of_threads)

makeFunctionalPrediction(path_to_otu_table = query_otu_table, path_to_reference_data = pwd_ref_data, path_to_temp_folder = pwd_op_folder, database_mode = "Ref99NR", normalize_by_copy_number = norm_by_cn, min_identity_to_reference = iden, normalize_pathways = norm_path)


## T2_D30_deseq2
# modify the following 8 lines as needed
query_otu_seq   = '~/Documents/R_analysis/jki_seq5/data_jki_seq5/Tax4Fun2/T2_D30_filt_deseq2.fasta'              # input OTU sequence file
query_otu_table = '~/Documents/R_analysis/jki_seq5/data_jki_seq5/Tax4Fun2/T2_D30_filt_deseq2.txt'          # input OTU table
pwd_op_folder   = '~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tax4Fun2_jki_seq5/T2_D30_deseq2'     # output directory
pwd_ref_data    = '~/Documents/R_analysis/jki_seq5/data_jki_seq5/Tax4Fun2/Tax4Fun2_ReferenceData_v2/'  # path to the default database (i.e., Tax4Fun2_ReferenceData_v2), need to be decompressed before use
norm_by_cn      = TRUE                          # normalize_by_copy_number (TRUE or FALSE)
norm_path       = TRUE                          # normalize_pathways (TRUE or FALSE)
iden            = 0.97                          # min_identity_to_reference, please modify as needed
num_of_threads  = 1                             # number of CPU cores to use, please modify as needed

# predict functions
runRefBlast(path_to_otus = query_otu_seq, path_to_reference_data = pwd_ref_data, path_to_temp_folder = pwd_op_folder, database_mode = "Ref99NR", use_force = T, num_threads = num_of_threads)

makeFunctionalPrediction(path_to_otu_table = query_otu_table, path_to_reference_data = pwd_ref_data, path_to_temp_folder = pwd_op_folder, database_mode = "Ref99NR", normalize_by_copy_number = norm_by_cn, min_identity_to_reference = iden, normalize_pathways = norm_path)

## T2_D60_deseq2
# modify the following 8 lines as needed
query_otu_seq   = '~/Documents/R_analysis/jki_seq5/data_jki_seq5/Tax4Fun2/T2_D60_filt_deseq2.fasta'              # input OTU sequence file
query_otu_table = '~/Documents/R_analysis/jki_seq5/data_jki_seq5/Tax4Fun2/T2_D60_filt_deseq2.txt'          # input OTU table
pwd_op_folder   = '~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tax4Fun2_jki_seq5/T2_D60_deseq2'     # output directory
pwd_ref_data    = '~/Documents/R_analysis/jki_seq5/data_jki_seq5/Tax4Fun2/Tax4Fun2_ReferenceData_v2/'  # path to the default database (i.e., Tax4Fun2_ReferenceData_v2), need to be decompressed before use
norm_by_cn      = TRUE                          # normalize_by_copy_number (TRUE or FALSE)
norm_path       = TRUE                          # normalize_pathways (TRUE or FALSE)
iden            = 0.97                          # min_identity_to_reference, please modify as needed
num_of_threads  = 1                             # number of CPU cores to use, please modify as needed

# predict functions
runRefBlast(path_to_otus = query_otu_seq, path_to_reference_data = pwd_ref_data, path_to_temp_folder = pwd_op_folder, database_mode = "Ref99NR", use_force = T, num_threads = num_of_threads)

makeFunctionalPrediction(path_to_otu_table = query_otu_table, path_to_reference_data = pwd_ref_data, path_to_temp_folder = pwd_op_folder, database_mode = "Ref99NR", normalize_by_copy_number = norm_by_cn, min_identity_to_reference = iden, normalize_pathways = norm_path)


## The end! : )







############
# ##Making functional predictions using the own database
# 
# #Generating the own reference database
# 
# library(Tax4Fun2)
# 
# # modify the following 4 lines as needed
# pwd_ref_data   = 'Tax4Fun2_ReferenceData_v2'   # path to Tax4Fun2's default database, need to be decompressed before use
# pwd_user_data  = 'genome_folder'               # path to the folder that holds the reference genome/MAG files
# name_user_data = 'name_of_user_database'       # specify the name of the generated database, specify only the name, do not include path here!
# gnm_ext        = 'fna'                         # extension of the genome files
# 
# # Generate your own database
# extractSSU(genome_folder = pwd_user_data, file_extension = gnm_ext, path_to_reference_data = pwd_ref_data)
# assignFunction(genome_folder = pwd_user_data, file_extension = gnm_ext, path_to_reference_data = pwd_ref_data, num_of_threads = num_threads, fast = TRUE)
# generateUserData(path_to_reference_data = pwd_ref_data, path_to_user_data = pwd_user_data, name_of_user_data = name_user_data, SSU_file_extension = "_16SrRNA.ffn", KEGG_file_extension = "_funPro.txt")
# 
# 
# 
# # Perform functional predictions
# library(Tax4Fun2)
# 
# # modify the following 10 lines as needed
# query_otu_seq   = "~/Documents/R_analysis/jki_seq5/data_jki_seq5/Tax4Fun2/T1_D30_filt.fasta"              # input OTU sequence file
# query_otu_table = "~/Documents/R_analysis/jki_seq5/data_jki_seq5/Tax4Fun2/T2_D30_filt.csv"          # input OTU table
# pwd_op_folder   = "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tax4Fun2_jki_seq5"      # output directory
# pwd_ref_data    = "~/Documents/R_analysis/jki_seq5/data_jki_seq5/Tax4Fun2/Tax4Fun2_ReferenceData_v2/"   # need to be the same as in step 1
# pwd_user_data   = "~/Documents/R_analysis/jki_seq5/data_jki_seq5/Tax4Fun2/Tax4Fun2_ReferenceData_v2/genome_folder"               # need to be the same as in step 1
# name_user_data  = 'Tax4Fun2_ReferenceData_v2'       # need to be the same as in step 1, specify only the name, do not include path here
# norm_by_cn      = TRUE                          # normalize_by_copy_number (TRUE or FALSE)
# norm_path       = TRUE                          # normalize_pathways (TRUE or FALSE)
# iden            = 0.97                          # min_identity_to_reference, please modify as needed
# num_of_threads  = 6                             # number of CPU cores to use , please modify as needed
# 
# # predict functions
# runRefBlast(path_to_otus = query_otu_seq, path_to_reference_data = pwd_ref_data, path_to_temp_folder = pwd_op_folder, database_mode = "Ref99NR", use_force = T, num_threads = num_of_threads, include_user_data = T, path_to_user_data = pwd_user_data, name_of_user_data = name_user_data)
# makeFunctionalPrediction(path_to_otu_table = query_otu_table, path_to_reference_data = pwd_ref_data, path_to_temp_folder = pwd_op_folder, database_mode = "Ref99NR", normalize_by_copy_number = norm_by_cn, min_identity_to_reference = iden, normalize_pathways = norm_path, include_user_data = T, path_to_user_data = pwd_user_data, name_of_user_data = name_user_data)