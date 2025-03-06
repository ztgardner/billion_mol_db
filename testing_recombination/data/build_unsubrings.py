import pandas as pd
import hashlib
import selfies as sf

def hash_smiles(smiles):
    """Generates a SHA-256 hash of a SMILES string."""
    return hashlib.sha256(smiles.encode()).hexdigest()

def process_ring_db(parquet_file, output_file):
    """Creates unsub_rings.parquet from ring_db.parquet with the required modifications."""
    try:
        # Load the original ring_db.parquet
        df = pd.read_parquet(parquet_file)

        # Create new DataFrame with the required columns
        new_df = pd.DataFrame()
        new_df["smiles"] = df["smiles"]  # Keep the SMILES column from ring_db
        new_df["id"] = new_df["smiles"].apply(hash_smiles)  # Generate hash for the ID
        new_df["id_ring"] = df["id"]  # Keep the original ID as id_ring
        new_df["id_sidechain"] = [[]] * len(new_df)  # Set id_sidechain as an empty list
        new_df["selfies"] = new_df["smiles"].apply(lambda x: sf.encoder(x) if pd.notna(x) else None)  # Convert SMILES to SELFIES

        # Remove duplicate rows based on the 'id' column
        new_df = new_df.drop_duplicates(subset=["id"])

        # Save the modified DataFrame as unsub_rings.parquet
        new_df.to_parquet(output_file, index=False)

        print(f"✅ Successfully created {output_file} with required modifications.")

    except Exception as e:
        print(f"❌ Error processing {parquet_file}: {e}")

if __name__ == "__main__":
    process_ring_db("ring_db.parquet", "unsub_rings.parquet")
