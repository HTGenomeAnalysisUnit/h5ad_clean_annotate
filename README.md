# h5ad clean and annotate

A small tool that help in cleaning, subsetting and annotation h5ad file efficiently.

This relies on bioalpha package from Bioturing for efficient processing, but no GPU is required.

## Order of operations

1. Clean index
2. Make new cell I
3. Add barcode based annotations
4. Add column based annotations
5. Rename columns using the rename map provided
6. Filter obs columns based on include/exclude lists
7. Sanitize obs columns
8. Set include_bc column based on the subset_bc list
