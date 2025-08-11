# h5ad clean and annotate

A small tool that help in cleaning, subsetting and annotation h5ad file efficiently.

The tool relies on bioalpha package from Bioturing for efficient processing, but no GPU is required.

## Usage

```bash
usage: clean_and_annotate.py [-h] --h5ad H5AD --out OUT --config CONFIG

Clean, annotate, subset and other useful operations on h5ad files.

options:
  -h, --help       show this help message and exit
  --h5ad H5AD      h5ad input file
  --out OUT        Output h5ad file
  --config CONFIG  JSON file defining configuration for cleaning and processing
```

## JSON config file

The operations to perform can be defined in a JSON file. A template structure is below and in the the `config_template.json` file.

```json5
{
	// Template to construct new cell IDs from obs columns. 
	// Any obs column can be accessed using curly brakets, e.g. {mycolumn}.
	"new_cell_id": "{tranche.id}--{tranche.name}--{index}", 
	
	// If true, the index will be cleaned by removing any --[0-9]+$ suffix.
	"clean_index": true,
	
	// A list of obs columns to retain
	"select_obs_columns": [
		"col1",
		"col2"
	],
	
	// A list of obs columns to remove
	"exclude_obs_columns": [
		"col1",
		"col2"
	],

	// If true, the obs column names will be sanitized by removing any special characters and spaces. 
	// >=  and <= will be replaced with greaterthan and lessthan strings
	"sanitize_obs_column_names": true,
	
	// Path to a text file containing barcodes to include in the output h5ad file.
	"subset_bc": "subset_barcodes.txt",
	
	// Path to a TSV file with 2 or more columns: barcode and annotations.
	// Annotation columns will be added to obs as new columns using the barcode as the key.
	"annot_bc": "cell_annotations.tsv",
	
	// A layer to be used as the default X in the output h5ad file.
	"X_layer": "layer_name",
	"keys": [
		"uns",
		"obsm",
		"raw/X",
		"layers/layer1"
	],

	// A dicitonary of "old_name": "new_name" pairs to rename obs columns.
	"rename_columns": {
		"old_name1": "new_name1",
		"old_name2": "new_name2"
	},

	// A list of dictionaries defining TSV files that can be used to annotate samples.
	// The tool will annotate obs with the table_annotation_columns from the TSV file
	// Merging keys in the input table and obs are defined by table_key_column and obs_key_column, repectively.
	"sample_annotations": [
		{
			"filename": "sample_annotations.tsv",
			"table_key_column": "sample_id",
			"obs_key_column": "sample_id",
			"table_annotation_columns": [
				"tissue",
				"treatment"
			],
			"annotation_name": "myanno1"
		},
		{
			"filename": "another_annotation.tsv",
			"table_key_column": "sample_id",
			"obs_key_column": "sample_id",
			"table_annotation_columns": [
				"ancestry"
			],
			"annotation_name": "myanno2"
		}
	]
}
```

## Order of operations

1. Clean index
2. Make new cell ID
3. Add barcode-based annotations
4. Rename columns using the rename map provided
5. Create new columns based on existing ones as defined in `make_columns`
6. Add column-based annotations
7. Filter obs columns based on include/exclude lists
9. Sanitize obs columns
10. Filter barcodes based on the subset_bc list
