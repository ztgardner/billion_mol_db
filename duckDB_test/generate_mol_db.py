import duckdb
import numpy as np
import pandas as pd
import pyarrow as pa
from pathlib import Path
import pyarrow.parquet as pq
from multiprocessing import Pool, cpu_count


def new_smiles(*fragments):
    """Function to generate new SMILES string (dummy example, replace with actual logic)"""
    # TODO Replace this with actual logic to combine two fragments
    return ".".join(fragments)


def process_chunk(chunk):
    """Function to process a chunk of combinations"""
    # Create a new list to store the processed results
    processed_data = []

    # Iterate over combinations and generate the new SMILES
    for index, row in chunk.iterrows():
        fragments = [row[f'smiles_{i}'] for i in range(1, len(row) // 2 + 1)]
        ids = [row[f'_id_{i}'] for i in range(1, len(row) // 2 + 1)]
        new_smiles_str = new_smiles(*fragments)

        # Append to the processed list
        processed_data.append({
            'smiles': new_smiles_str,
            **{f'_id_{i}': ids[i - 1] for i in range(1, len(ids) + 1)}
        })

    return pd.DataFrame(processed_data)


# Function to parallelize the processing of data
def parallel_process(data: pd.DataFrame, func, n_cores=None):
    if n_cores is None:
        n_cores = cpu_count()
    print(f"Running {n_cores} cores for parallel data generation")
    data_split = np.array_split(data, n_cores)
    pool = Pool(n_cores)
    data = pd.concat(pool.map(func, data_split))
    pool.close()
    pool.join()
    return data


def main():
    # Variable and file paths
    HOME_DIR = Path(__file__).resolve().parent.parent
    fragment_file = HOME_DIR / 'data' / 'fragment_db.parquet'
    mol_db_file = HOME_DIR / 'data' / 'mol_db.parquet'
    id_name = "_id"
    smiles_name = "d1_smiles"

    # Load database fragments into DuckDB
    con = duckdb.connect(database=':memory:')
    con.execute(f"CREATE TABLE fragments AS SELECT * FROM read_parquet('{fragment_file}')")
    print(con.execute(f"SELECT * FROM fragments LIMIT 5").fetch_df())

    # Generate all combinations of 2 fragments, ensuring no duplicates (order does not matter)
    combinations_query_2 = f"""
    SELECT a.{id_name} as _id_1, a.{smiles_name} as smiles_1, 
           b.{id_name} as _id_2, b.{smiles_name} as smiles_2
    FROM fragments AS a, fragments AS b
    WHERE a.{id_name} <= b.{id_name}
    """
    # # Generate all combinations of up to 4 fragments, ensuring no duplicates (order does not matter)
    # combinations_query_4 = f"""
    # WITH combinations_2 AS (
    #     SELECT a.{id_name} as _id_1, a.{smiles_name} as smiles_1, 
    #            b.{id_name} as _id_2, b.{smiles_name} as smiles_2
    #     FROM fragments a, fragments b
    #     WHERE a.{id_name} < b.{id_name}
    # ),
    # combinations_3 AS (
    #     SELECT a.{id_name} as _id_1, a.{smiles_name} as smiles_1, 
    #            b.{id_name} as _id_2, b.{smiles_name} as smiles_2,
    #            c.{id_name} as _id_3, c.{smiles_name} as smiles_3
    #     FROM fragments a, fragments b, fragments c
    #     WHERE a.{id_name} < b.{id_name} AND b.{id_name} < c.{id_name}
    # ),
    # combinations_4 AS (
    #     SELECT a.{id_name} as _id_1, a.{smiles_name} as smiles_1, 
    #            b.{id_name} as _id_2, b.{smiles_name} as smiles_2,
    #            c.{id_name} as _id_3, c.{smiles_name} as smiles_3,
    #            d.{id_name} as _id_4, d.{smiles_name} as smiles_4
    #     FROM fragments a, fragments b, fragments c, fragments d
    #     WHERE a.{id_name} < b.{id_name} AND b.{id_name} < c.{id_name} AND c.{id_name} < d.{id_name}
    # )
    # SELECT * FROM combinations_2
    # UNION ALL
    # SELECT * FROM combinations_3
    # UNION ALL
    # SELECT * FROM combinations_4
    # """

    # Execute the query and fetch combinations in chunks to avoid memory issues
    combinations_df = con.execute(combinations_query_2).fetch_df()

    # Process the combinations in parallel
    processed_df = parallel_process(combinations_df, process_chunk)

    # Convert the processed data into a parquet file
    table = pa.Table.from_pandas(processed_df)
    pq.write_table(table, mol_db_file)

    print(f"Mol DB written to {mol_db_file} with {processed_df.shape[0]} rows")


if __name__ == "__main__":
    main()
