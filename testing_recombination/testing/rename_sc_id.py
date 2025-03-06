# import pandas as pd
# import string
#
#
# def generate_labels(n):
#     """Generate sequential labels A, B, ..., Z, AA, AB, ... for n items."""
#     labels = []
#     alphabet = list(string.ascii_uppercase)
#
#     for i in range(n):
#         label = ""
#         temp = i
#         while temp >= 0:
#             label = alphabet[temp % 26] + label
#             temp = temp // 26 - 1
#         labels.append(label)
#
#     return labels
#
#
# def rename_sidechain_ids(parquet_file):
#     """Renames the 'id' column in sidechain.parquet with sequential letters (A, B, C, ...)."""
#
#     # Load Parquet file
#     try:
#         df = pd.read_parquet(parquet_file)
#     except Exception as e:
#         print(f"❌ Error loading {parquet_file}: {e}")
#         return
#
#     # Generate new labels
#     num_rows = len(df)
#     new_labels = generate_labels(num_rows)
#
#     # Assign new labels to 'id' column
#     df["id"] = new_labels
#
#     # Save back to Parquet
#     try:
#         df.to_parquet(parquet_file, index=False)
#         print(f"✅ Successfully renamed 'id' column in {parquet_file}")
#     except Exception as e:
#         print(f"❌ Error saving {parquet_file}: {e}")
#
#
# if __name__ == "__main__":
#     parquet_path = "sidechain.parquet"  # Adjust path if needed
#     rename_sidechain_ids(parquet_path)

import pandas as pd

def remove_rows_from_parquet(parquet_file):
    """Removes the 25th and 27th rows from the given Parquet file and saves the modified version."""
    try:
        # Load the Parquet file
        df = pd.read_parquet(parquet_file)

        # Drop rows with indices 24 (25th row) and 26 (27th row)
        df = df.drop(index=[24, 26], errors='ignore')

        # Save the modified DataFrame back to the same Parquet file
        df.to_parquet(parquet_file, index=False)

        print("✅ Successfully removed the 25th and 27th rows.")

    except Exception as e:
        print(f"❌ Error processing {parquet_file}: {e}")

if __name__ == "__main__":
    remove_rows_from_parquet("sidechain.parquet")
