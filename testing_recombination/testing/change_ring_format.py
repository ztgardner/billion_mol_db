import pandas as pd
import re

def format_nodes(value):
    """Formats node strings: (index, element) → (index, 'element')"""
    if isinstance(value, str):
        formatted = re.sub(r"\((\d+),\s*([A-Za-z\*]+)\)", r"(\1, '\2')", value)
        return formatted
    return value

def format_edges(value):
    """Formats edge strings: (n1, n2, bond) → (n1, n2, 'bond')"""
    if isinstance(value, str):
        formatted = re.sub(r"\((\d+),\s*(\d+),\s*([A-Za-z]+)\)", r"(\1, \2, '\3')", value)
        return formatted
    return value

def fix_ring_format(input_file, output_file):
    """Reads the parquet file, formats nodes and edges, and saves the corrected version."""
    print(f"Loading {input_file}...")
    df = pd.read_parquet(input_file)

    print("Fixing node and edge formatting...")

    df["nodes"] = df["nodes"].astype(str).apply(format_nodes)
    df["edges"] = df["edges"].astype(str).apply(format_edges)

    print(f"Saving fixed file to {output_file}...")
    df.to_parquet(output_file, index=False)
    print("✅ File saved successfully!")

if __name__ == "__main__":
    fix_ring_format("../data/ring_db.parquet", "ring_db.parquet")
