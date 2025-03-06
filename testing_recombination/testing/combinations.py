import sys
import duckdb
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from pathlib import Path
import itertools
import os

class GenComs_Ring_Sidechains:
    def __init__(self, procnumber, output_dir="output"):
        """
        Initializes the class with processor-specific ring and sidechain datasets.
        """
        self.procnumber = procnumber
        self.ring_file = f"ring_{procnumber}.parquet"  # Processor-specific ring file
        self.sidechain_file = "sidechain.parquet"  # Shared sidechain file
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.output_file = self.output_dir / f"combinations_{procnumber}.parquet"
        self.con = duckdb.connect()

    def load_data(self):
        """
        Loads ring and sidechain datasets into DuckDB tables.
        """
        if not Path(self.ring_file).exists():
            raise FileNotFoundError(f"Processor {self.procnumber}: {self.ring_file} not found.")
        if not Path(self.sidechain_file).exists():
            raise FileNotFoundError(f"Processor {self.procnumber}: {self.sidechain_file} not found.")

        print(f"Processor {self.procnumber}: Loading {self.ring_file} and {self.sidechain_file} into DuckDB...")

        # Load into DuckDB
        self.con.execute(f"CREATE TABLE rings AS SELECT * FROM '{self.ring_file}'")
        self.con.execute(f"CREATE TABLE sidechains AS SELECT * FROM '{self.sidechain_file}'")

    def ensure_valid_parquet(self):
        """
        Ensures that the output Parquet file is valid by creating an empty file with the correct schema if needed.
        """
        parquet_schema = pa.schema([
            ("id", pa.int32()),  # Unique combination ID
            ("id_ring", pa.string()),  # Ring ID
            ("id_sidechain", pa.list_(pa.string())),  # List of sidechain IDs
            ("smiles", pa.string()),  # Placeholder
            ("selfies", pa.string()),  # Placeholder
        ])

        if not self.output_file.exists() or os.path.getsize(self.output_file) == 0:
            print(f"Processor {self.procnumber}: Initializing {self.output_file} with an empty valid schema.")
            empty_table = pa.Table.from_pandas(pd.DataFrame(columns=parquet_schema.names), schema=parquet_schema)
            pq.write_table(empty_table, self.output_file)

    def generate_combinations(self):
        """
        Processes rings one at a time, generating all possible sidechain combinations.
        """
        self.ensure_valid_parquet()

        # Define schema inside this function to avoid NameError
        parquet_schema = pa.schema([
            ("id", pa.int32()),  # Unique combination ID
            ("id_ring", pa.string()),  # Ring ID
            ("id_sidechain", pa.list_(pa.string())),  # List of sidechain IDs (fixed)
            ("smiles", pa.string()),  # Placeholder
            ("selfies", pa.string()),  # Placeholder
        ])

        # Read existing combinations to determine the starting ID
        try:
            existing_table = pq.read_table(self.output_file)
            existing_df = existing_table.to_pandas()
            next_id = existing_df["id"].max() + 1 if not existing_df.empty else 0
        except:
            next_id = 0

        # Read rings and sidechains tables into DuckDB
        rings_df = self.con.execute("SELECT * FROM rings").fetchdf()
        sidechains_df = self.con.execute("SELECT * FROM sidechains").fetchdf()

        print(f"Processor {self.procnumber}: Columns in rings_df: {rings_df.columns.tolist()}")
        print(f"Processor {self.procnumber}: Columns in sidechains_df: {sidechains_df.columns.tolist()}")

        required_ring_columns = ["id", "sub_points"]
        required_sidechain_columns = ["id"]

        for col in required_ring_columns:
            if col not in rings_df.columns:
                raise KeyError(f"Processor {self.procnumber}: Column '{col}' not found in rings_df.")

        for col in required_sidechain_columns:
            if col not in sidechains_df.columns:
                raise KeyError(f"Processor {self.procnumber}: Column '{col}' not found in sidechains_df.")

        # Ensure sidechain IDs are stored as proper **strings**
        sidechains_df["id"] = sidechains_df["id"].astype(str)

        # Process one ring at a time
        for _, ring in rings_df.iterrows():
            ring_id = ring["id"]
            n_subs = ring["sub_points"]

            # Generate all possible combinations of `n_subs` sidechains (including repetitions)
            sidechain_combinations = list(itertools.product(sidechains_df.itertuples(index=False), repeat=n_subs))

            results = []
            for sidechain_set in sidechain_combinations:
                # Fix: Ensure `sub_ids` are stored as **strings**, not ASCII lists
                sub_ids = [str(sc.id) for sc in sidechain_set]  # <-- FIXED

                row = {
                    "id": next_id,  # Unique combination ID
                    "id_ring": str(ring_id),
                    "id_sidechain": sub_ids,  # List of **string** sidechain IDs
                    "smiles": "N/A",  # Placeholder
                    "selfies": "N/A",
                }
                results.append(row)
                next_id += 1  # Increment the ID

            df_combinations = pd.DataFrame(results, columns=[
                "id", "id_ring", "id_sidechain", "smiles", "selfies",
            ])

            if not df_combinations.empty:
                existing_table = pq.read_table(self.output_file)
                new_table = pa.Table.from_pandas(df_combinations, schema=parquet_schema)
                combined_table = pa.concat_tables([existing_table, new_table])
                pq.write_table(combined_table, self.output_file)

            print(f"Processor {self.procnumber}: Processed ring {ring_id} with {n_subs} substituent positions.")

    def close_connection(self):
        """Closes the DuckDB connection."""
        self.con.close()


# Run script with processor-specific input
if __name__ == "__main__":
    try:
        procnumber = int(sys.argv[1])  # Get processor number from CLI argument
    except (IndexError, ValueError):
        print("Error: Missing or invalid processor number argument. Exiting.")
        sys.exit(1)

    output_dir = "output"  # Store files in an output directory

    # Initialize and run the generator
    generator = GenComs_Ring_Sidechains(procnumber, output_dir)
    generator.load_data()
    generator.generate_combinations()
    generator.close_connection()
