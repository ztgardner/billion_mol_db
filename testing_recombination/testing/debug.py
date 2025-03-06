# import duckdb
#
# con = duckdb.connect()
# sidechains_df = con.execute("SELECT * FROM 'sidechain.parquet'").fetchdf()
# print(sidechains_df.head())
# print(sidechains_df.dtypes)  # Check if 'id' is being treated as an object, string, or something else


# import pandas as pd
#
# df = pd.read_parquet("output/combinations_1.parquet")
# print(df.head())  # See how IDs are stored
# print(df.dtypes)  # Verify column types



import duckdb

rings_file = "ring_1.parquet"  # Modify based on your test
con = duckdb.connect()

# Fetch first row of rings.parquet
ring_data = con.execute(f"SELECT id, nodes, edges FROM read_parquet('{rings_file}') LIMIT 5").fetchall()
con.close()

print("Sample rows from rings.parquet:")
for row in ring_data:
    print("ID:", row[0])
    print("Nodes:", row[1])  # Expected: [(1, 'C'), (2, 'N')]
    print("Edges:", row[2])  # Expected: [(1, 2, 'Single'), (2, 3, 'Double')]
    print("=" * 50)
