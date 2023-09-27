



# ---- get file ----
if (!hasArg(file)) file <- file.choose()

file <- "Data/..."


# ---- import fasta file ----
fasta <- Biostrings::readAAStringSet(filepath = file)




