import argparse
import bioalpha as bsc
import pandas as pd
from string import Formatter
import json

# Example of sample_annotations JSON file:
# [
#   {
#	 "filename": "sample_annotations.tsv",
#	 "table_key_column": "sample_id",
#	 "obs_key_column": "sample_id",
#	 "table_annotation_columns": ["tissue", "treatment"], # Columns to add to obs
#	 "annotation_name": "myanno1" # This define a suffix that is added to table_column name when a column with the same name is alredy in obs. If omitted annotation is added as suffix
#   },
#   {
#	 "filename": "another_annotation.tsv",
#	 "table_key_column": "sample_id",
#	 "obs_key_column": "sample_id",
#	 "table_annotation_columns": ["ancestry"],
#	 "annotation_name": "myanno2"
#   }
# ]

def main():
	#Cmd line args
	parser = argparse.ArgumentParser(description='Clean cell IDs and eventually subset the dataset')
	parser.add_argument('--h5ad', action='store', required=True,
						help='h5ad input file')
	parser.add_argument('--out', action='store', required=True,
						help='Output h5ad file')
	parser.add_argument('--new_cell_id', action='store', default=None, 
					 	help="A pattern like '{col1}--{col2}--{index}' on how the new cell IDs should look like. If not provided, the old cell IDs are used.")
	parser.add_argument('--subset_bc', action='store', default=None,
						help='File with a list of barcodes to subset')
	parser.add_argument('--annot_bc', action='store', default=None,
						help='Tab-separated file with header. Cell barcodes are expected in column cell_id, other cols will be added to obs.')
	parser.add_argument('--annot_samples', action='store', default=None,
						help='JSON file defining sample annotations tables and strategy. If provided, the sample annotations will be added to obs.')
	parser.add_argument('--layer', action='store', default=None,
						help='A layer to use as X in the output. If not provided, the input X is used.')
	parser.add_argument('--keys', action='store', default=None,
						help='Comma-separated list of features to keep in new file. Ex: “uns,obsm”. If not provided, no features are saved to output.')
	parser.add_argument( '--rename_col', type=str, action='append',
						help='Specify a column to rename in "old_name,new_name" format. Can be used multiple times.')
	parser.add_argument('--clean_index', action='store_true',
						help='If provided remove --* suffix from original index before using it in new_cell_id pattern.')
	args = parser.parse_args()

	# Load the h5ad file
	adata = bsc.h5ad_map.H5ADMap(args.h5ad)
	print(f'Loaded {args.h5ad} with {adata.n_obs} cells and {adata.n_vars} genes.')

	# Initialize an empty dictionary to store the renaming map
	rename_map = {}

	# Check if the --rename_col argument was provided
	if args.rename_col:
		print(f"Processing {len(args.rename_col)} column rename pairs...")
		# Iterate over the list of "old_name,new_name" strings
		for pair in args.rename_col:
			# Split the string by the comma
			parts = pair.split(',')
			# Ensure there are exactly two parts
			if len(parts) != 2:
				raise ValueError(f"Invalid format. Expected 'old_name,new_name', but got '{pair}'")
			
			# Trim whitespace from both parts
			old_name = parts[0].strip()
			new_name = parts[1].strip()

			# Check for empty names
			if not old_name or not new_name:
				raise ValueError(f"Column names cannot be empty in pair '{pair}'")

			# Add the pair to our dictionary
			rename_map[old_name] = new_name
			print(f"  - Mapping '{old_name}' -> '{new_name}'")

	# If a layer is specified, use it as X
	outlayer = 'X'
	if args.layer is not None:
		if args.layer not in adata.layers:
			raise ValueError(f"Layer '{args.layer}' not found in the input file. Available layers: {list(adata.layers.keys())}")
		print(f'Using layer "{args.layer}" as X.')
		outlayer = f'layers/{args.layer}'

	# If keys are specified, make a list, otherwise set to None
	keys = None
	if args.keys is not None:
		keys = args.keys.split(',')
		print(f'Keeping the following features in output: {keys}')
	else:
		print('WARN - No additonal feature key (--keys) specified. No uns, obsm, etc. will be saved to output.')

	# If clean_index is specified, remove the --* suffix from the index
	if args.clean_index:
		print("Cleaning index by removing --* suffix")
		adata.obs.index = adata.obs.index.str.replace(r'--\d+$', '', regex=True)

	# Given a pattern string like "{col1}--{col2}--{index}" create a new cell ID accordingly
	if args.new_cell_id is not None:
		# Build the new cell IDs vectorially. This is much faster than apply() and
		# works with dots in column names (e.g. 'tranche.id').
		print(f'Create new cell IDs with pattern: {args.new_cell_id}')
		new_ids = pd.Series([""] * adata.n_obs, index=adata.obs.index, dtype=str)
		for literal_text, field_name, _, _ in Formatter().parse(args.new_cell_id):
			if literal_text:
				new_ids += literal_text
			if field_name:
				if field_name == 'index':
					values = adata.obs.index
				else:
					values = adata.obs[field_name]
				new_ids += values.astype(str)
		adata.obs['cell_id'] = new_ids.values
		adata.obs.set_index('cell_id', inplace=True)
		print(f'Example of new IDs: {adata.obs.index[:5].tolist()}')
		print(f'Cell ID updated')
	
	# Set a include_bc flag in obs if the cell ID is in the subset_bc file
	adata.obs['include_bc'] = True  # Default to include all cells
	if args.subset_bc is not None:
		print(f'Reading subset barcodes from file: {args.subset_bc}')
		subset_bc = pd.read_csv(args.subset_bc, header=None, names=['cell_id'])['cell_id'].astype(str).tolist()
		adata.obs['include_bc'] = adata.obs.index.isin(subset_bc)
		print(f'N barcodes present in subset: {adata.obs["include_bc"].sum()} out of {adata.n_obs} total cells.')
	
	# If annot_bc is provided, read it and merge with adata.obs
	if args.annot_bc is not None:
		print(f'Reading annotation file: {args.annot_bc}')
		annot_bc = pd.read_csv(args.annot_bc, sep='\t')
		if 'cell_id' not in annot_bc.columns:
			raise ValueError("The annotation file must contain a 'cell_id' column.")
		annot_bc.set_index('cell_id', inplace=True)
		# Check if there is at least some overlap between adata.obs.index and annot_bc['cell_id']
		if not adata.obs.index.isin(annot_bc.index).any():
			raise ValueError("No matching cell IDs found between adata.obs.index and annot_bc['cell_id'].")
		adata.obs = adata.obs.join(annot_bc, how='left')
		print(f'Annotation file merged. New obs columns: {list(adata.obs.columns)}')

	# If annot_samples is provided, read the files defined in the JSON and merge with adata.obs based on configured columns
	if args.annot_samples is not None:
		print(f'Reading annotation samples config: {args.annot_samples}')
		with open(args.annot_samples, 'r') as f:
			annot_samples = json.load(f)
		for annotation_config in annot_samples:
			print(f'processing configuration: {annotation_config}')
			# annotation config is a dist with filename, annotation_name, table_column, obs_column
			filename = annotation_config['filename']
			table_key_column = annotation_config['table_key_column']
			obs_key_column = annotation_config['obs_key_column']
			annotation_columns = annotation_config['table_annotation_columns']
			annotation_name = annotation_config.get('annotation_name', 'annotation')
			
			print(f'Reading annotation file: {filename}')
			annot_samples_df = pd.read_csv(filename, sep='\t')

			full_table_columns = annotation_columns + [table_key_column]

			if not set(full_table_columns).issubset(annot_samples_df.columns):
				raise ValueError(f"One of the defined columns is not found in the annotation file '{filename}'.")
			if obs_key_column not in adata.obs.columns:
				raise ValueError(f"Column '{obs_key_column}' not found in adata.obs.")
			# Merge the annotation file with adata.obs based on the specified columns
			print(f'Merging annotation file {filename} with adata.obs on {obs_key_column} and {table_key_column}')
			right_df = annot_samples_df[full_table_columns].set_index(table_key_column)
			adata.obs = adata.obs.join(right_df, 
										on=obs_key_column, how='left', 
										rsuffix=f'_{annotation_name}')
			
			# If table_key_column is now in obs, remove it
			if table_key_column in adata.obs.columns:
				adata.obs.drop(columns=[table_key_column], inplace=True)

			# Check if new added column contains NaN values since this is not compatible with anndata
			# In this case column type is set to string, and NaN is replaced with 'NOT_ASSIGNED'
			for col in annotation_columns:
				if adata.obs[col].isnull().any():
					print(f"\nColumn '{col}' in obs contains NaN values. These will be replaced with 'NOT_ASSIGNED'.")
					adata.obs[col] = adata.obs[col].astype(str).fillna('NOT_ASSIGNED')
			
			print(f'Annotation file {filename} merged. New obs columns: {list(adata.obs.columns)}')
	
	# Rename columns in obs if rename_map is provided
	if len(rename_map) > 0:
		print(f'Renaming columns in obs according to the provided map: {rename_map}')
		for old_name, new_name in rename_map.items():
			if old_name in adata.obs.columns:
				adata.obs.rename(columns={old_name: new_name}, inplace=True)
			else:
				print(f"Warning: Column '{old_name}' not found in obs. Skipping renaming.")

	# Diet subset using include_bc flag
	print("Subset file")
	if keys is None:
		adata = adata.diet_subset(args.out, X=outlayer, row_index="include_bc")
	else:
		adata = adata.diet_subset(args.out, X=outlayer, row_index="include_bc", keys=keys)

	print("All done")

if __name__ == '__main__':
	main()
