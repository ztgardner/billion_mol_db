import duckdb
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from multiprocessing import Pool, cpu_count

fragment_file = 'fragment_db.parquet'
mol_db_file = 'mol_db.parquet'

# Load database A into DuckDB
con = duckdb.connect(database=':memory:')
con.execute(f"CREATE TABLE fragments AS SELECT * FROM read_parquet('{fragment_file}')")


# Function to generate new SMILES string (dummy example, replace with actual logic)
def new_smiles(fragment_a, fragment_b):
    # TODO Replace this with actual logic to combine two fragments
    return fragment_a + "." + fragment_b


# Function to process a chunk of combinations
def process_chunk(chunk):
    # Create a new list to store the processed results
    processed_data = []

    # Iterate over combinations and generate the new SMILES
    for index, row in chunk.iterrows():
        fragment_a = row['smiles_a']
        fragment_b = row['smiles_b']
        _id_a = row['_id_a']
        _id_b = row['_id_b']
        new_smiles_str = new_smiles(fragment_a, fragment_b)

        # Append to the processed list
        processed_data.append({
            'smiles': new_smiles_str,
            '_id_a': _id_a,
            '_id_b': _id_b
        })

    return pd.DataFrame(processed_data)


# Function to parallelize the processing of data
def parallel_process(data, func, n_cores=None):
    if n_cores is None:
        n_cores = cpu_count()
    data_split = np.array_split(data, n_cores)
    pool = Pool(n_cores)
    data = pd.concat(pool.map(func, data_split))
    pool.close()
    pool.join()
    return data


# Generate all combinations of fragments
combinations_query = """
SELECT a._id as _id_a, a.smiles as smiles_a, 
       b._id as _id_b, b.smiles as smiles_b 
FROM fragments a, fragments b
"""

# Execute the query and fetch combinations in chunks to avoid memory issues
con.execute(combinations_query)
combinations_df = con.fetch_df()

# Process the combinations in parallel
processed_df = parallel_process(combinations_df, process_chunk)

# Convert the processed data into a parquet file
table = pa.Table.from_pandas(processed_df)
pq.write_table(table, mol_db_file)

print(f"Mol DB written to {mol_db_file}")
