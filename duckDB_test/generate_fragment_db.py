import pandas as pd
import pyarrow as pa
from pathlib import Path
import pyarrow.parquet as pq
HOME_DIR = Path(__file__).resolve().parent.parent

# Input and output file paths
csv_file = HOME_DIR / 'data' / 'data_10mil_actual.csv'  # Path to your CSV file
parquet_file = HOME_DIR / 'data' / 'fragment_db.parquet'  # Path where the output Parquet file will be saved

# Read the CSV into a Pandas DataFrame
df = pd.read_csv(csv_file)

# Convert DataFrame to an Arrow Table
table = pa.Table.from_pandas(df)

# Write the Arrow Table to a Parquet file
pq.write_table(table, parquet_file)

print(f"Parquet file '{parquet_file}' created successfully for dataset of {df.shape[0]} rows!")
