# h5ad clean and annotate

A small tool that help in cleaning, subsetting and annotation h5ad file efficiently.

This relies on bioalpha package from Bioturing for efficient processing, but no GPU is required.

## Order of operations

1. Clean index
2. Make new cell ID
3. Filter barcodes based on subset_bc list
4. Add barcode based annotations
5. Add column based annotations
6. Rename columns using the rename map provided
7. Filter obs columns based on include/exclude lists
8. Sanitize obs columns
