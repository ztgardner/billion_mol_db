import duckdb
import pandas as pd
import concurrent.futures
import pyarrow as pa
from pathlib import Path
import pyarrow.parquet as pq
import os
import shutil


def new_smiles(fragment_a, fragment_b):
    """Logic to combine two SMILES strings."""
    return fragment_a + "." + fragment_b  # Example combination logic


# Function to process and write each batch to its own file
def process_and_write_batch(df_batch, parquet_schema, batch_id, temp_dir):
    df_batch["smiles"] = df_batch.apply(lambda row: new_smiles(row["smiles_a"], row["smiles_b"]), axis=1)
    df_batch = df_batch[["id_a", "id_b", "smiles"]]  # Select only the relevant columns

    # Convert DataFrame to Apache Arrow table
    table = pa.Table.from_pandas(df_batch, schema=parquet_schema)

    # Write to a temporary Parquet file
    temp_file = os.path.join(temp_dir, f"batch_{batch_id}.parquet")
    pq.write_table(table, temp_file)


def main():
    # Variable and file paths
    HOME_DIR = Path(__file__).resolve().parent.parent
    fragment_file = HOME_DIR / 'data' / 'fragment_db.parquet'
    mol_db_file = HOME_DIR / 'data' / 'mol_db.parquet'
    temp_dir = HOME_DIR / 'temp_batches'  # Directory to store temporary batch files
    temp_dir.mkdir(exist_ok=True)

    id_name = "_id"
    smiles_name = "d1_smiles"
    batch_size = 1000  # Adjust this batch size as necessary

    # Load the molecule fragments from a file into DuckDB
    con = duckdb.connect()
    con.execute(f"CREATE TABLE fragments AS SELECT * FROM '{fragment_file}'")  # Load fragments into DuckDB

    # Set up Parquet schema
    parquet_schema = pa.schema([("id_a", pa.string()),
                                ("id_b", pa.string()),
                                ("smiles", pa.string())])

    # Create the query for unique combinations
    query = f"""
        SELECT a.{id_name} as id_a, a.{smiles_name} as smiles_a, 
               b.{id_name} as id_b, b.{smiles_name} as smiles_b
        FROM fragments AS a, fragments AS b
        WHERE a.{id_name} <= b.{id_name}
        """

    # Execute the query and use fetch_df_chunk to fetch data in chunks
    cursor = con.execute(query)

    executor = concurrent.futures.ProcessPoolExecutor()
    batch_id = 0

    while True:
        # Fetch the next chunk of data using fetch_df_chunk()
        df_batch = cursor.fetch_df_chunk(batch_size)
        if df_batch.empty:
            break  # No more data

        # Process each batch in parallel and write to its own file
        future = executor.submit(process_and_write_batch, df_batch=df_batch, parquet_schema=parquet_schema,
                                 batch_id=batch_id, temp_dir=temp_dir)

        batch_id += 1

    # After processing, merge all temp files into the final file
    with pq.ParquetWriter(mol_db_file, parquet_schema) as writer:
        for temp_file in temp_dir.glob("*.parquet"):
            table = pq.read_table(temp_file)
            writer.write_table(table)

    total_rows = con.execute(f"SELECT COUNT(*) FROM ({query})").fetchone()[0]
    print(f"Mol DB with {total_rows} rows written to {mol_db_file}.")

    # Cleanup temporary files and connections
    shutil.rmtree(temp_dir)  # Remove entire temp directory
    con.close()
    executor.shutdown()


if __name__ == "__main__":
    main()
