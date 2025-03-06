import pandas as pd
import glob
import datetime


def append_parquet_files():
    """Appends all combinations_proc.parquet files into one big parquet file with a timestamp."""

    # Generate timestamp in format YYYY_M_D_HHMM
    timestamp = datetime.datetime.now().strftime("%Y_%m_%d_%H%M")
    output_file = f"combined_{timestamp}.parquet"

    # Find all combinations_proc.parquet files in the directory
    parquet_files = glob.glob("output/combinations_*.parquet")

    if not parquet_files:
        print("❌ No combinations_proc.parquet files found.")
        return

    # Read and append all Parquet files
    dataframes = []
    for file in parquet_files:
        try:
            df = pd.read_parquet(file)
            dataframes.append(df)
            print(f"✅ Successfully loaded {file}")
        except Exception as e:
            print(f"❌ Error loading {file}: {e}")

    # Concatenate all dataframes
    if dataframes:
        combined_df = pd.concat(dataframes, ignore_index=True)

        # Save to new Parquet file
        try:
            combined_df.to_parquet(output_file, index=False)
            print(f"🚀 Successfully saved combined file: {output_file}")
        except Exception as e:
            print(f"❌ Error saving combined file: {e}")
    else:
        print("⚠️ No valid data to combine.")


if __name__ == "__main__":
    append_parquet_files()
