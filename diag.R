# diag.R

cat("Activating renv...\n")
source("renv/activate.R")

# Print library paths
sink("diag_libpaths.txt")
cat("LIBPATHS IN BATCH:\n")
print(.libPaths())
sink()

cat("Done.\n")

