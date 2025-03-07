import sys
import duckdb
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import ast
from rdkit import Chem
import selfies as sf
import time
import hashlib

def update_smiles_hash_in_parquet(procnumber):
    """Computes a hash of the SMILES column and updates the 'id' column in combinations_proc.parquet."""

    combinations_file = f"output/combinations_{procnumber}.parquet"

    # Load Parquet file into Pandas
    try:
        df = pd.read_parquet(combinations_file)
    except Exception as e:
        print(f"❌ Error loading {combinations_file}: {e}")
        return

    # Ensure 'smiles' column exists
    if "smiles" not in df.columns:
        print(f"❌ 'smiles' column not found in {combinations_file}. Skipping update.")
        return

    # Compute SHA-256 hash of each SMILES string and update the 'id' column
    def hash_smiles(smiles):
        return hashlib.sha256(smiles.encode()).hexdigest() if isinstance(smiles, str) else None

    df["id"] = df["smiles"].apply(hash_smiles)

    # Write the modified DataFrame back to Parquet
    try:
        df.to_parquet(combinations_file, index=False)
        print(f"✅ SMILES hashes updated successfully in 'id' column for {combinations_file}")
    except Exception as e:
        print(f"❌ Error writing back to {combinations_file}: {e}")

def remove_duplicate_rows(procnumber):
    """Removes duplicate rows based on the 'id' column in combinations_proc.parquet."""

    combinations_file = f"output/combinations_{procnumber}.parquet"

    # Load Parquet file into Pandas
    try:
        df = pd.read_parquet(combinations_file)
    except Exception as e:
        print(f"❌ Error loading {combinations_file}: {e}")
        return

    # Remove duplicates based on 'id' column while keeping the first occurrence
    df = df.drop_duplicates(subset=["id"], keep="first")

    # Write the cleaned DataFrame back to Parquet
    try:
        df.to_parquet(combinations_file, index=False)
        print(f"✅ Duplicate rows removed successfully from {combinations_file}")
    except Exception as e:
        print(f"❌ Error writing back to {combinations_file}: {e}")

def count_rows(procnumber):
    """Counts the number of rows in combinations_{procnumber}.parquet using DuckDB."""
    combinations_file = f"output/combinations_{procnumber}.parquet"
    con = duckdb.connect()

    try:
        result = con.execute(f"SELECT COUNT(*) FROM read_parquet('{combinations_file}')").fetchone()
        con.close()
        return result[0] if result else 0
    except Exception as e:
        print(f"❌ Error counting rows in {combinations_file}: {e}")
        con.close()
        return 0

def process_n_rows(procnumber):
    """Processes all rows from combinations_{procnumber}.parquet and updates their SMILES and SELFIES."""
    start_time = time.time()  # Start the timer
    combinations_file = f"output/combinations_{procnumber}.parquet"

    # Get the number of rows dynamically
    n = count_rows(procnumber)
    if n == 0:
        print(f"⚠️ No rows to process in {combinations_file}. Exiting.")
        return

    con = duckdb.connect()

    # Load the data
    try:
        df = con.execute(f"SELECT * FROM read_parquet('{combinations_file}')").df()
    except Exception as e:
        print(f"❌ Error loading {combinations_file}: {e}")
        return

    # Process each row
    for row_index in range(n):
        row_start_time = time.time()  # Start time for this row

        ring_id = df.iloc[row_index]["id_ring"]
        sidechain_ids = df.iloc[row_index]["id_sidechain"]

        print(f"\n🔄 Processing row {row_index + 1}/{n}:")
        print(f"   🧩 Ring ID: {ring_id}")
        print(f"   🛠️ Sidechain IDs: {sidechain_ids}")

        # Load ring and sidechain data
        ring_df = pd.read_parquet(f"ring_{procnumber}.parquet")
        sidechain_df = pd.read_parquet("sidechain.parquet")

        # Extract graphs **with the correct row index**
        ring_graph, sidechain_graphs, _, _ = extract_data(procnumber, ring_df, sidechain_df, row_index)

        if not ring_graph or not sidechain_graphs:
            print(f"❌ Skipping row {row_index + 1}, missing valid graph data.")
            continue

        # Combine into final graph
        final_graph = combine_graphs(ring_graph, sidechain_graphs)
        if not final_graph:
            print(f"❌ Skipping row {row_index + 1}, failed to construct final graph.")
            continue

        # Convert graph to SMILES
        smiles_string = graph_to_smiles(final_graph)
        if smiles_string:
            print(f"✅ Final SMILES: {smiles_string}")
            update_smiles_in_parquet(procnumber, smiles_string, row_index)

            # Convert to SELFIES
            selfies_string = smiles_to_selfies(smiles_string)
            update_selfies_in_parquet(procnumber, selfies_string, row_index)
        else:
            print(f"⚠️ SMILES conversion failed for row {row_index + 1}, skipping update.")

        # Row processing time
        row_end_time = time.time()
        print(f"⏱️ Row {row_index + 1} processed in {row_end_time - row_start_time:.2f} seconds.")

    con.close()

    # Total processing time
    end_time = time.time()
    print(f"\n🚀 Finished processing {n} rows in {end_time - start_time:.2f} seconds.")


def load_combination(procnumber):
    """Loads and prints the first row from the combinations file for a given processor number."""
    combinations_file = f"output/combinations_{procnumber}.parquet"

    # Connect to DuckDB
    con = duckdb.connect()

    # Fetch a single row
    try:
        combination_data = con.execute(f"SELECT * FROM read_parquet('{combinations_file}') LIMIT 1").fetchone()
    except Exception as e:
        print(f"Error reading {combinations_file}: {e}")
        return None

    if not combination_data:
        print(f"No data found in {combinations_file}.")
        return None

    # Print the raw row for debugging
    print(f"Processor {procnumber} - Combination Row:")
    print(combination_data)

    con.close()
    return combination_data

def build_graph(nodes, edges):
    """Builds a NetworkX graph from node and edge lists."""
    G = nx.Graph()

    # Add nodes with labels
    for node_idx, label in nodes:
        G.add_node(node_idx, label=label)

    # Add edges with bond labels
    for node1, node2, bond in edges:
        G.add_edge(node1, node2, bond=bond)

    # Print graph details
    print(f"Graph Nodes: {list(G.nodes(data=True))}")
    print(f"Graph Edges: {list(G.edges(data=True))}")

    return G

def visualize_graph(G, title):
    """Visualizes a NetworkX graph with node labels."""
    plt.figure(figsize=(5, 5))
    pos = nx.spring_layout(G)  # Positioning of nodes
    labels = nx.get_node_attributes(G, 'label')
    edge_labels = {(u, v): d["bond"] for u, v, d in G.edges(data=True)}

    nx.draw(G, pos, with_labels=True, labels=labels, node_color="lightblue", edge_color="gray")
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_color="red")
    plt.title(title)
    plt.show()

def safe_parse_nodes_edges(value, label="unknown"):
    """Safely parses node or edge lists from a string representation, handling errors."""
    if isinstance(value, str):
        try:
            parsed_value = ast.literal_eval(value)  # Use safer parsing
            if isinstance(parsed_value, list):
                return parsed_value  # Ensure it’s a list before returning
        except (SyntaxError, ValueError):
            print(f"Error: Could not parse {label} - {value}")  # Debug message

    return []  # Return empty list on failure

def extract_data(procnumber, ring_df, sidechain_df, row_index):
    """
    Extracts the ring and sidechains for the specified row index.
    Returns:
      - ring_graph: NetworkX graph of the ring
      - sidechain_graphs: List of NetworkX graphs of sidechains
      - ring_id: The ID of the ring
      - sidechain_ids: List of IDs of the sidechains
    """

    # Load the full combinations file (since we're working with row indices now)
    combinations_file = f"output/combinations_{procnumber}.parquet"
    con = duckdb.connect()

    try:
        df = con.execute(f"SELECT * FROM read_parquet('{combinations_file}')").df()
    except Exception as e:
        print(f"❌ Error loading {combinations_file}: {e}")
        return None, None, None, None

    # Make sure row index is valid
    if row_index >= len(df):
        print(f"❌ Invalid row index {row_index}, skipping...")
        return None, None, None, None

    # Get the correct row
    ring_id = df.iloc[row_index]["id_ring"]
    sidechain_ids = df.iloc[row_index]["id_sidechain"]

    print(f"🔎 Processing row {row_index + 1}: Ring ID = {ring_id}, Sidechain IDs = {sidechain_ids}")

    # Ensure sidechain_ids is a list
    if isinstance(sidechain_ids, str):
        try:
            sidechain_ids = ast.literal_eval(sidechain_ids)  # Convert from string representation
        except (SyntaxError, ValueError):
            print(f"⚠️ Warning: Could not parse sidechain IDs: {sidechain_ids}")
            sidechain_ids = []

    # Ensure ring_df["id"] is treated as a string
    ring_df["id"] = ring_df["id"].astype(str)
    sidechain_df["id"] = sidechain_df["id"].astype(str)

    # Extract the ring row
    ring_row = ring_df.loc[ring_df["id"] == ring_id]
    if ring_row.empty:
        print(f"⚠️ Ring {ring_id} NOT found in ring_{procnumber}.parquet! Skipping...")
        return None, None, ring_id, sidechain_ids

    print(f"✅ Ring Found: {ring_row.to_dict(orient='records')}")

    # Parse ring nodes and edges safely
    ring_nodes = safe_parse_nodes_edges(ring_row.iloc[0]["nodes"], label="ring nodes")
    ring_edges = safe_parse_nodes_edges(ring_row.iloc[0]["edges"], label="ring edges")

    # Build the ring graph
    ring_graph = build_graph(ring_nodes, ring_edges)

    # Extract and build sidechain graphs
    sidechain_graphs = []
    for i, sid in enumerate(sidechain_ids):
        sidechain_row = sidechain_df.loc[sidechain_df["id"] == sid]

        if sidechain_row.empty:
            print(f"⚠️ Sidechain {i + 1} (ID: {sid}) NOT found in sidechain.parquet! Skipping.")
            continue

        print(f"✅ Sidechain {i + 1}: {sidechain_row.to_dict(orient='records')}")

        sidechain_nodes = safe_parse_nodes_edges(sidechain_row.iloc[0]["nodes"], label=f"sidechain {i+1} nodes")
        sidechain_edges = safe_parse_nodes_edges(sidechain_row.iloc[0]["edges"], label=f"sidechain {i+1} edges")

        sidechain_graph = build_graph(sidechain_nodes, sidechain_edges)
        sidechain_graphs.append(sidechain_graph)

    return ring_graph, sidechain_graphs, ring_id, sidechain_ids

def combine_graphs(ring_graph, sidechain_graphs):
    """Combines the ring graph and all sidechains, connecting the first '*' in each."""

    # Create a new combined graph
    combined_graph = nx.Graph()

    # Track node indices for merging
    node_offset = 0  # Offset to prevent duplicate indices when merging graphs
    node_map = {}  # To track renumbered nodes

    # Add the ring graph first
    for node, data in ring_graph.nodes(data=True):
        new_idx = node + node_offset
        combined_graph.add_node(new_idx, label=data["label"])
        node_map[node] = new_idx

    for n1, n2, bond in ring_graph.edges(data=True):
        combined_graph.add_edge(node_map[n1], node_map[n2], bond=bond["bond"])

    # Find first starred atom in the ring
    ring_star_nodes = [n for n, d in combined_graph.nodes(data=True) if '*' in d['label']]

    if not ring_star_nodes:
        print("❌ No star nodes found in the ring graph! Aborting combination.")
        return None

    # Track the index shift
    node_offset = max(combined_graph.nodes()) + 1

    # Add each sidechain
    for sidechain_graph in sidechain_graphs:
        sidechain_star_nodes = [n for n, d in sidechain_graph.nodes(data=True) if '*' in d['label']]

        if not sidechain_star_nodes:
            print("⚠️ No star nodes found in a sidechain, skipping this sidechain.")
            continue

        # Get the first starred atom in both the ring and sidechain
        ring_star = ring_star_nodes.pop(0)
        sidechain_star = sidechain_star_nodes[0]

        # Offset the sidechain node indices
        sidechain_node_map = {}
        for node, data in sidechain_graph.nodes(data=True):
            new_idx = node + node_offset
            combined_graph.add_node(new_idx, label=data["label"])
            sidechain_node_map[node] = new_idx

        for n1, n2, bond in sidechain_graph.edges(data=True):
            combined_graph.add_edge(sidechain_node_map[n1], sidechain_node_map[n2], bond=bond["bond"])

        # Add the bond between the first stars in ring and sidechain
        combined_graph.add_edge(ring_star, sidechain_node_map[sidechain_star], bond="SINGLE")

        # Remove a single '*' from labels
        combined_graph.nodes[ring_star]['label'] = combined_graph.nodes[ring_star]['label'].replace('*', '', 1)
        combined_graph.nodes[sidechain_node_map[sidechain_star]]['label'] = \
        combined_graph.nodes[sidechain_node_map[sidechain_star]]['label'].replace('*', '', 1)

        # Update the offset
        node_offset = max(combined_graph.nodes()) + 1

    return combined_graph

def draw_graph(G, title="Final Combined Graph"):
    """Draws a NetworkX graph with labeled nodes and edges."""
    plt.figure(figsize=(8, 8))
    pos = nx.spring_layout(G)  # Positioning of nodes
    labels = nx.get_node_attributes(G, 'label')
    edge_labels = {(u, v): d["bond"] for u, v, d in G.edges(data=True)}

    nx.draw(G, pos, with_labels=True, labels=labels, node_color="lightblue", edge_color="gray")
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_color="red")
    plt.title(title)
    plt.show()

def graph_to_smiles(G):
    """Converts a NetworkX graph into an RDKit SMILES string with error handling."""
    try:
        mol = Chem.RWMol()
        atom_map = {}

        # Step 1: Add atoms to RDKit molecule
        for idx, data in G.nodes(data=True):
            atom_label = data["label"].replace("*", "")  # Remove any remaining '*'
            atom = Chem.Atom(atom_label)
            atom_idx = mol.AddAtom(atom)
            atom_map[idx] = atom_idx

        # Step 2: Define bond types
        bond_dict = {
            "SINGLE": Chem.BondType.SINGLE,
            "DOUBLE": Chem.BondType.DOUBLE,
            "TRIPLE": Chem.BondType.TRIPLE,
            "AROMATIC": Chem.BondType.AROMATIC
        }

        # Step 3: Add bonds to RDKit molecule
        for n1, n2, bond_data in G.edges(data=True):
            bond_type = bond_dict.get(bond_data["bond"], Chem.BondType.SINGLE)  # Default to SINGLE
            mol.AddBond(atom_map[n1], atom_map[n2], bond_type)

        # Step 4: Generate and return SMILES
        smiles = Chem.MolToSmiles(mol, canonical=True)
        return smiles

    except Exception as e:
        print(f"❌ Conversion to SMILES failed: {e}")
        return None  # Return None instead of crashing

def smiles_to_selfies_and_update(procnumber, ring_id, sidechain_ids):
    """Converts SMILES to SELFIES and updates the correct row in combinations_proc.parquet."""

    combinations_file = f"output/combinations_{procnumber}.parquet"
    con = duckdb.connect()

    sidechain_str = str(sidechain_ids)

    # Debug: Print what we're searching for
    print(f"🔎 Searching for row: ring_id = {ring_id}, sidechains = {sidechain_str}")

    # Fetch the corresponding SMILES
    try:
        smiles_row = con.execute(
            f"SELECT smiles FROM read_parquet('{combinations_file}') WHERE id_ring = ? AND id_sidechain = ?",
            [ring_id, sidechain_str]
        ).fetchone()

        if not smiles_row:
            print(
                f"⚠️ No matching row found for ring {ring_id} and sidechains {sidechain_str}. Skipping SELFIES update.")
            return

        smiles = smiles_row[0]
    except Exception as e:
        print(f"❌ Error retrieving SMILES: {e}")
        return

    # Convert SMILES to SELFIES
    try:
        selfies_str = sf.encoder(smiles)
    except Exception as e:
        print(f"❌ SELFIES conversion failed for SMILES {smiles}: {e}")
        selfies_str = None

    if not selfies_str:
        print(f"⚠️ SELFIES conversion failed. Skipping update.")
        return

    # Update the SELFIES column directly in DuckDB
    try:
        con.execute(
            f"""
            UPDATE read_parquet('{combinations_file}')
            SET selfies = ?
            WHERE id_ring = ? AND id_sidechain = ?
            """,
            [selfies_str, ring_id, sidechain_str]
        )
        print(f"✅ SELFIES updated successfully in {combinations_file}")
    except Exception as e:
        print(f"❌ Error updating SELFIES in {combinations_file}: {e}")

    con.close()

def update_smiles_in_parquet(procnumber, new_smiles, row_index):
    """Updates the SMILES column for a specific row index in combinations_proc.parquet using Pandas."""

    combinations_file = f"output/combinations_{procnumber}.parquet"

    # Load Parquet file into Pandas
    try:
        df = pd.read_parquet(combinations_file)
    except Exception as e:
        print(f"❌ Error loading {combinations_file}: {e}")
        return

    # Ensure the row index is within bounds
    if row_index >= len(df):
        print(f"⚠️ Row index {row_index} is out of bounds. Skipping update.")
        return

    # Update the SMILES column for the given row index
    df.at[row_index, "smiles"] = new_smiles

    # Write the modified DataFrame back to Parquet
    try:
        df.to_parquet(combinations_file, index=False)
        print(f"✅ SMILES updated successfully for row {row_index + 1}")
    except Exception as e:
        print(f"❌ Error writing back to {combinations_file}: {e}")

def update_selfies_in_parquet(procnumber, new_selfies, row_index):
    """Updates the SELFIES column for a specific row index in combinations_proc.parquet using Pandas."""

    combinations_file = f"output/combinations_{procnumber}.parquet"

    # Load Parquet file into Pandas
    try:
        df = pd.read_parquet(combinations_file)
    except Exception as e:
        print(f"❌ Error loading {combinations_file}: {e}")
        return

    # Ensure the row index is within bounds
    if row_index >= len(df):
        print(f"⚠️ Row index {row_index} is out of bounds. Skipping update.")
        return

    # Update the SELFIES column for the given row index
    df.at[row_index, "selfies"] = new_selfies

    # Write the modified DataFrame back to Parquet
    try:
        df.to_parquet(combinations_file, index=False)
        print(f"✅ SELFIES updated successfully for row {row_index + 1}")
    except Exception as e:
        print(f"❌ Error writing back to {combinations_file}: {e}")

def smiles_to_selfies(smiles):
    """Converts SMILES to SELFIES with error handling."""
    try:
        selfies_str = sf.encoder(smiles)
        return selfies_str
    except Exception as e:
        print(f"❌ SELFIES conversion failed: {e}")
        return None  # Return None if conversion fails

if __name__ == "__main__":
    try:
        procnumber = int(sys.argv[1])
    except (IndexError, ValueError):
        print("Error: Invalid arguments. Usage: python builder.py <procnumber> <num_rows>")
        sys.exit(1)

    process_n_rows(procnumber)
    update_smiles_hash_in_parquet(procnumber)
    remove_duplicate_rows(procnumber)