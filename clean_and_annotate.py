import argparse
import bioalpha as bsc
import pandas as pd
from string import Formatter
import json

# Example of sample_annotations JSON file:
# {
# "new_cell_id": "{tranche.id}--{tranche.name}--{index}",
# "clean_index": true, # If true, remove --* suffix from original index before using it in new_cell_id pattern.
# "select_obs_columns": ['col1','col2'], # List of columns to select from obs
# "exclude_obs_columns": ['col1','col2'], # List of columns to remove from obs, ignored if select_obs_columns is provided
# "sanitize_obs_column_names": true # If true, sanitize obs column names
# "subset_bc": "subset_barcodes.txt", # File with a list of barcodes to subset
# "annot_bc": "cell_annotations.tsv", # Tab-separated file with header containing cell annotations. Cell barcodes are expected in column cell_id, other cols will be added to obs.
# "X_layer": "layer_name", # A layer to use as X in the output. If not provided, the input X is used.
# "keys": ["uns","obsm", "raw/X", "layers/layer1"], # List of additional features to keep in output. If not provided, no features are saved to output.
# "rename_columns": {
# 	"old_name1": "new_name1"
# 	"old_name2": "new_name2"
# }, # Specify columns to rename as "old_name": "new_name"
# "sample_annotations": [
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
# }

VERSION = "0.1.0"

def main():
	#Cmd line args
	parser = argparse.ArgumentParser(description='Clean cell IDs and eventually subset the dataset')
	parser.add_argument('--h5ad', action='store', required=True,
						help='h5ad input file')
	parser.add_argument('--out', action='store', required=True,
						help='Output h5ad file')
	parser.add_argument('--config', action='store', default=None, required=True,
						help='JSON file defining configuration for cleaning and processing')
	args = parser.parse_args()

	print(f"== H5AD CLEAN AND ANNOTATE VERSION {VERSION} ==")

	# Load the configuration file
	try:
		with open(args.config, 'r') as f:
			config = json.load(f)
	except FileNotFoundError:
		raise FileNotFoundError(f"Configuration file '{args.config}' not found.")
	except json.JSONDecodeError:
		raise ValueError(f"Configuration file '{args.config}' is not a valid JSON file.")

	# Log the configuration
	print(f'Configuration loaded from {args.config}:')
	for scope_name, value in config.items():
		print(f'  {scope_name}: {value}')

	print("== PREPARING TO PROCESS ==")

	# Set rename_map dictionary to config['rename_columns'] or an empty dict if does not exist
	if not isinstance(config['rename_columns'], dict):
		raise ValueError("The 'rename_columns' configuration must be a dictionary with 'old_name': 'new_name' pairs.")
	rename_map = config.get('rename_columns', {})
	print(f"- Found {len(rename_map)} columns to rename")

	# If a layer is not specified, just use X
	outlayer = config.get('X_layer', 'X')
	print(f'- Using layer "{outlayer}" as X in output.')

	# If keys are specified, make a list, otherwise set to None
	keys = config.get('keys', None)
	if keys is not None:
		print(f'- Keeping the following features in output: {keys}')
	else:
		print('WARN - No additonal feature keys specified. No uns, obsm, layers will be saved to output.')

	# If subset_bc is specified read the file and store the barcodes
	subset_bc_file = config.get('subset_bc', None)
	subset_bc = pd.DataFrame()
	if subset_bc_file is not None:
		print(f'- Reading subset barcodes from file: {subset_bc_file}')
		subset_bc = pd.read_csv(subset_bc_file, header=None, names=['cell_id'])['cell_id'].astype(str).tolist()

	# If annot_bc is specified read the file and store the barcode annotations
	annot_bc_file = config.get('annot_bc', None)
	annot_bc = pd.DataFrame()
	if annot_bc_file is not None:
		print(f'- Reading annotation file: {annot_bc_file}')
		annot_bc = pd.read_csv(annot_bc_file, sep='\t')
		if 'cell_id' not in annot_bc.columns:
			raise ValueError("The cell annotation file must contain a 'cell_id' column.")
		annot_bc.set_index('cell_id', inplace=True)

	sanitize_col_names = config.get('sanitize_obs_column_names', False)
	print(f'- Sanitize column names: {sanitize_col_names}')
	clean_index = config.get('clean_index', False)
	print(f'- Clean index: {clean_index}')
	new_cell_id = config.get('new_cell_id', None)
	if new_cell_id is not None: print(f'- New cell ID pattern: {new_cell_id}')

	# Prepare list of columns to select/exclude
	select_obs_columns = config.get('select_obs_columns', [])
	exclude_obs_columns = config.get('exclude_obs_columns', [])
	if len(select_obs_columns) > 0:
		print(f'- Selecting the following columns from obs: {select_obs_columns}')
		if len(exclude_obs_columns) > 0:
			print("Warning: 'select_obs_columns' and 'exclude_obs_columns' are both specified. 'select_obs_columns' will be used.")
	elif len(exclude_obs_columns) > 0:
		print(f'- Excluding the following columns from obs: {exclude_obs_columns}')

	print ("== START PROCESSING ==")
	# Load the h5ad file
	adata = bsc.h5ad_map.H5ADMap(args.h5ad)
	print(f'Loaded {args.h5ad} with {adata.n_obs} cells and {adata.n_vars} genes.')

	# If a layer is specified, use it as X
	if outlayer != 'X' and outlayer not in adata.layers:
		raise ValueError(f"Layer '{outlayer}' not found in the input file. Available layers: {list(adata.layers.keys())}")
	print(f'Using layer "{outlayer}" as X.')
	outlayer = f'layers/{outlayer}'

	# If clean_index is specified, remove the --* suffix from the index
	if clean_index:
		print("Cleaning index by removing --* suffix")
		adata.obs.index = adata.obs.index.str.replace(r'--\d+$', '', regex=True)

	# Given a pattern string like "{col1}--{col2}--{index}" create a new cell ID accordingly
	if new_cell_id is not None:
		# Build the new cell IDs vectorially. This is much faster than apply() and
		# works with dots in column names (e.g. 'tranche.id').
		print(f'Create new cell IDs with pattern: {new_cell_id}')
		new_ids = pd.Series([""] * adata.n_obs, index=adata.obs.index, dtype=str)
		for literal_text, field_name, _, _ in Formatter().parse(new_cell_id):
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
	if not subset_bc.empty:
		adata.obs['include_bc'] = adata.obs.index.isin(subset_bc)
		print(f'N barcodes present in subset: {adata.obs["include_bc"].sum()} out of {adata.n_obs} total cells.')
	
	# If annot_bc is provided, read it and merge with adata.obs
	if not annot_bc.empty:
		print("Annotating cells with provided cell annotations.")
		# Check if there is at least some overlap between adata.obs.index and annot_bc['cell_id']
		if not adata.obs.index.isin(annot_bc.index).any():
			raise ValueError("No matching cell IDs found between adata.obs.index and annot_bc['cell_id'].")
		adata.obs = adata.obs.join(annot_bc, how='left')
		print(f'Annotation file merged. New obs columns: {list(adata.obs.columns)}')

	# If annot_samples is provided, read the files defined in the JSON and merge with adata.obs based on configured columns
	annot_samples = config.get('sample_annotations', [])
	if len(annot_samples) > 0:
		print(f'Found {len(annot_samples)} sample annotation configurations to process.')
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
		
		print(f'Annotation file {filename} merged')

	# Rename columns in obs if rename_map is provided
	if len(rename_map) > 0:
		print(f'Renaming columns in obs according to the provided map')
		for old_name, new_name in rename_map.items():
			if old_name in adata.obs.columns:
				adata.obs.rename(columns={old_name: new_name}, inplace=True)
			else:
				print(f"Warning: Column '{old_name}' not found in obs. Skipping renaming.")

	# Get the list of columns to keep in obs
	if len(exclude_obs_columns) > 0 or len(select_obs_columns) > 0:
		final_column_set = set(adata.obs.columns) - set(exclude_obs_columns)
		if len(select_obs_columns) > 0:
			final_column_set = set(select_obs_columns)
		print(f'Reducing obs columns to: {final_column_set}')
		
		# Check all columns in final_column_set are in adata.obs
		if not final_column_set.issubset(adata.obs.columns):
			raise ValueError(f"Some columns to select are not found in adata.obs: {set(select_obs_columns) - set(adata.obs.columns)}")

		# Select the columns from obs
		adata.obs = adata.obs[list(final_column_set)].copy()

	# If sanitize_col_names is True, sanitize obs column names
	if sanitize_col_names:
		original_col_names = set(adata.obs.columns)
		print("Sanitizing obs column names")
		# Replace space, dot, column, semicolon, slash, and dash with underscore
		adata.obs.columns = adata.obs.columns.str.replace(r'[ .;,:-/]', '_', regex=True)
		# Remove leading and trailing underscores/spaces
		adata.obs.columns = adata.obs.columns.str.strip('_')
		adata.obs.columns = adata.obs.columns.str.strip(' ')
		# Replace >= and > with greaterthan
		adata.obs.columns = adata.obs.columns.str.replace(r'>=', 'greaterthan', regex=True)
		adata.obs.columns = adata.obs.columns.str.replace(r'>', 'greaterthan', regex=True)
		# Replave <= and < with lessthan
		adata.obs.columns = adata.obs.columns.str.replace(r'<=', 'lessthan', regex=True)
		adata.obs.columns = adata.obs.columns.str.replace(r'<', 'lessthan', regex=True)
		# Remove brackets and curly brackets
		adata.obs.columns = adata.obs.columns.str.replace(r'[\[\]{}]', '', regex=True)
		# Set everything to lowercase
		adata.obs.columns = adata.obs.columns.str.lower()

		# Print number of col names modified
		modified_col_names = set(adata.obs.columns) - original_col_names
		print(f"Sanitized {len(modified_col_names)}")
		print(f"Modified column names: {', '.join(modified_col_names)}")

	print("== FINISHED PROCESSING ==")

	# Diet subset using include_bc flag
	print("Saving processed file")
	if keys is None:
		adata = adata.diet_subset(args.out, X=outlayer, row_index="include_bc")
	else:
		adata = adata.diet_subset(args.out, X=outlayer, row_index="include_bc", keys=keys)

	print("== ALL DONE ==")

if __name__ == '__main__':
	main()
