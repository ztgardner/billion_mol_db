import duckdb
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from pathlib import Path
import itertools
import os


class GenComs_Ring_Sidechains:
    def __init__(self, ring_file, sidechain_file, output_dir):
        """
        Initializes the class with input ring and sidechain datasets and the output directory.
        The Parquet file name will be automatically determined to avoid overwrites.
        """
        self.ring_file = ring_file
        self.sidechain_file = sidechain_file
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.output_file = self.get_unique_filename()
        self.con = duckdb.connect()

    def get_unique_filename(self):
        """
        Generate a unique Parquet filename by appending _1, _2, etc., if a file already exists.
        """
        base_filename = "combinations"
        extension = ".parquet"
        output_path = self.output_dir / f"{base_filename}{extension}"

        counter = 1
        while output_path.exists():
            output_path = self.output_dir / f"{base_filename}_{counter}{extension}"
            counter += 1

        return output_path

    def load_data(self):
        """
        Loads ring and sidechain datasets into DuckDB tables.
        """
        self.con.execute(f"CREATE TABLE rings AS SELECT * FROM '{self.ring_file}'")
        self.con.execute(f"CREATE TABLE sidechains AS SELECT * FROM '{self.sidechain_file}'")

    def ensure_valid_parquet(self):
        """
        Ensures that the output Parquet file is valid by creating an empty file with the correct schema if needed.
        """
        parquet_schema = pa.schema([
            ("ring_id", pa.string()),
            ("ring_smiles", pa.string()),
            ("sub_id", pa.list_(pa.string())),  # Ordered list of sidechain IDs
            ("sub_smiles", pa.list_(pa.string())),  # Ordered list of sidechain SMILES
            ("combination", pa.string())  # Placeholder for now
        ])

        # If the file does not exist or is empty, create a valid empty Parquet file
        if not os.path.exists(self.output_file) or os.path.getsize(self.output_file) == 0:
            print(f"Initializing {self.output_file} with an empty valid schema.")
            empty_table = pa.Table.from_pandas(pd.DataFrame(columns=parquet_schema.names), schema=parquet_schema)
            pq.write_table(empty_table, self.output_file)

    def generate_combinations(self):
        """
        Processes rings one at a time, generating all possible sidechain combinations
        (including repeated sidechains) based on the number of substitution points.
        Saves results to the Parquet file iteratively.
        """
        # Ensure output file is valid before writing
        self.ensure_valid_parquet()

        # Read rings and sidechains tables into DuckDB
        rings_df = self.con.execute("SELECT * FROM rings").fetchdf()
        sidechains_df = self.con.execute("SELECT * FROM sidechains").fetchdf()

        # Define the output Parquet schema
        parquet_schema = pa.schema([
            ("ring_id", pa.string()),
            ("ring_smiles", pa.string()),
            ("sub_id", pa.list_(pa.string())),  # Ordered list of sidechain IDs
            ("sub_smiles", pa.list_(pa.string())),  # Ordered list of sidechain SMILES
            ("combination", pa.string())  # Placeholder for now
        ])

        # Process one ring at a time
        for _, ring in rings_df.iterrows():
            ring_id = ring["id"]
            ring_smiles = ring["smiles"]
            n_subs = ring["n_subs"]

            # Generate all possible combinations of `n_subs` sidechains (including repetitions)
            sidechain_combinations = list(itertools.product(sidechains_df.itertuples(index=False), repeat=n_subs))

            results = []
            for sidechain_set in sidechain_combinations:
                sub_ids = [sc.id for sc in sidechain_set]
                sub_smiles = [sc.smiles for sc in sidechain_set]

                # Construct a row with ordered lists
                row = {
                    "ring_id": str(ring_id),
                    "ring_smiles": str(ring_smiles),
                    "sub_id": sub_ids,  # List of sidechain IDs
                    "sub_smiles": sub_smiles,  # List of sidechain SMILES
                    "combination": "N/A"  # Placeholder
                }
                results.append(row)

            # Convert results to a DataFrame
            df_combinations = pd.DataFrame(results,
                                           columns=["ring_id", "ring_smiles", "sub_id", "sub_smiles", "combination"])

            # Append to Parquet file manually
            if not df_combinations.empty:
                # Read the existing Parquet file
                existing_table = pq.read_table(self.output_file)

                # Convert DataFrame to Arrow Table with the same schema
                new_table = pa.Table.from_pandas(df_combinations, schema=parquet_schema)

                # Concatenate existing and new tables
                combined_table = pa.concat_tables([existing_table, new_table])

                # Write back to Parquet file
                pq.write_table(combined_table, self.output_file)

            print(f"Processed ring {ring_id} with {n_subs} substituent positions.")

    def close_connection(self):
        """
        Closes the DuckDB connection.
        """
        self.con.close()


# Example usage
if __name__ == "__main__":
    # Define file paths
    HOME_DIR = Path.cwd()
    ring_file = HOME_DIR / "data" / "ring_list.parquet"
    sidechain_file = HOME_DIR / "data" / "sidechain_list.parquet"
    output_dir = HOME_DIR / "output"  # Store files in an output directory

    # Initialize and run the class
    generator = GenComs_Ring_Sidechains(ring_file, sidechain_file, output_dir)
    generator.load_data()
    generator.generate_combinations()
    generator.close_connection()
