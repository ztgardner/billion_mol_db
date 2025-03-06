import pandas as pd


def modify_ring_parquet():
    """Loads ring.parquet, keeps only the first row, and sets sub_points to 1."""

    input_file = "ring.parquet"  # Adjust path if needed
    output_file = "ring_modified.parquet"  # Save modified version

    try:
        # Load the Parquet file
        df = pd.read_parquet(input_file)

        if df.empty:
            print("❌ Error: ring.parquet is empty.")
            return

        # Keep only the first row
        df = df.iloc[:1]

        # Set sub_points to 1
        if "sub_points" in df.columns:
            df.loc[:, "sub_points"] = 1
        else:
            print("⚠️ Warning: 'sub_points' column not found in ring.parquet.")

        # Save the modified file
        df.to_parquet(output_file, index=False)
        print(f"✅ Modified ring.parquet saved as {output_file}")

    except Exception as e:
        print(f"❌ Error processing ring.parquet: {e}")


if __name__ == "__main__":
    modify_ring_parquet()
