pip install rdkit

import networkx as nx
from rdkit import Chem
import sqlite3
import pandas as pd
from collections import Counter

def create_database():
    conn = sqlite3.connect('molecules.db')
    cursor = conn.cursor()

    cursor.execute('''
    CREATE TABLE IF NOT EXISTS input_smiles (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        smiles TEXT NOT NULL
    )
    ''')

    cursor.execute('''
    CREATE TABLE IF NOT EXISTS subgraphs (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        parent_id INTEGER,
        smiles TEXT NOT NULL UNIQUE,
        nodes TEXT,
        edges TEXT,
        atoms INTEGER,
        FOREIGN KEY(parent_id) REFERENCES input_smiles(id)
    )
    ''')

    conn.commit()
    conn.close()

def add_atom_column_if_needed(conn, atom_type):
    cursor = conn.cursor()
    cursor.execute(f"PRAGMA table_info(subgraphs)")
    columns = [info[1] for info in cursor.fetchall()]

    sanitized_atom_type = f'atom_{atom_type}'  # Prefix atom type to ensure it's safe as a column name

    if sanitized_atom_type not in columns:
        cursor.execute(f"ALTER TABLE subgraphs ADD COLUMN {sanitized_atom_type} INTEGER DEFAULT 0")
        conn.commit()

def smiles_to_graph(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Warning: SMILES string '{smiles}' could not be parsed.")
        return None
    g = nx.Graph()
    for atom in mol.GetAtoms():
        g.add_node(atom.GetIdx(), element=atom.GetSymbol())
    for bond in mol.GetBonds():
        g.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond_type=bond.GetBondType())
    return g

def graph_to_mol(graph):
    mol = Chem.RWMol()
    node_to_idx = {}
    for node, data in graph.nodes(data=True):
        a = Chem.Atom(data['element'])
        mol_idx = mol.AddAtom(a)
        node_to_idx[node] = mol_idx
    try:
        for start, end, data in graph.edges(data=True):
            mol.AddBond(node_to_idx[start], node_to_idx[end], data['bond_type'])
        Chem.SanitizeMol(mol)
    except (Chem.KekulizeException, Chem.AtomValenceException) as e:
        return None
    return mol

def reindex_subgraph(subgraph, cut_vertex):
    mapping = {cut_vertex: 0}
    new_index = 1
    for node in subgraph.nodes():
        if node != cut_vertex:
            mapping[node] = new_index
            new_index += 1
    return nx.relabel_nodes(subgraph, mapping)

def identify_and_cut_bridges(graph):
    subgraphs = []
    bridge_bond_types = {}
    unique_smiles = set()
    original_graph = graph.copy()
    bridges = list(nx.bridges(original_graph))

    for bridge in bridges:
        bond_type = graph.edges[bridge]['bond_type']
        bridge_bond_types[bridge] = bond_type

        # Remove the bridge to generate subgraphs
        graph.remove_edge(*bridge)
        components = [graph.subgraph(c).copy() for c in nx.connected_components(graph)]

        for component in components:
            mol = graph_to_mol(component)
            if mol:
                smiles = Chem.MolToSmiles(mol, canonical=True)
                if smiles not in unique_smiles:
                    unique_smiles.add(smiles)
                    cut_vertex = bridge[0] if bridge[0] in component.nodes() else bridge[1]
                    reindexed_subgraph = reindex_subgraph(component, cut_vertex)
                    subgraphs.append((reindexed_subgraph, smiles))  # Store the reindexed subgraph and its canonical SMILES

        # Restore the bridge
        graph.add_edge(*bridge, bond_type=bond_type)

    return subgraphs, bridge_bond_types

def format_nodes_and_edges(graph):
    nodes = [(node, data['element']) for node, data in graph.nodes(data=True)]
    edges = []
    for u, v, data in graph.edges(data=True):
        bond_type = data['bond_type']
        bond_str = f"{u},{v}: {bond_type.name.capitalize()}"
        edges.append(bond_str)
    return nodes, edges

def count_atoms_from_nodes(nodes_str):
    nodes = eval(nodes_str)  # Safely convert the string back to a list
    atom_counts = Counter([element for _, element in nodes])
    return atom_counts

def write_subgraphs_to_db(parent_id, subgraphs):
    conn = sqlite3.connect('molecules.db')
    cursor = conn.cursor()

    for graph, smiles in subgraphs:
        nodes, edges = format_nodes_and_edges(graph)
        nodes_str = str(nodes)
        edges_str = ', '.join(edges)
        atoms_count = len(nodes)
        atom_counts = count_atoms_from_nodes(nodes_str)

        # Add new columns dynamically for any new atom types encountered
        for atom_type in atom_counts.keys():
            add_atom_column_if_needed(conn, atom_type)

        # Get the current list of columns in the table
        cursor.execute(f"PRAGMA table_info(subgraphs)")
        columns = [info[1] for info in cursor.fetchall()]

        # Prepare data for all atom columns, defaulting to 0 if not present
        atom_columns = []
        for atom in columns[6:]:
            sanitized_atom = f'atom_{atom}'
            atom_columns.append(atom_counts.get(atom[5:], 0))  # Remove 'atom_' prefix when counting

        # Construct the SQL INSERT statement with dynamic columns
        subgraph_data = (parent_id, smiles, nodes_str, edges_str, atoms_count, *atom_columns)
        placeholders = ', '.join(['?'] * len(subgraph_data))
        cursor.execute(f'INSERT OR IGNORE INTO subgraphs (parent_id, smiles, nodes, edges, atoms, {", ".join(columns[6:])}) VALUES ({placeholders})', subgraph_data)

    conn.commit()
    conn.close()

def process_smiles(smiles, parent_id):
    graph = smiles_to_graph(smiles)
    if graph is None:
        return  # Skip invalid SMILES
    subgraphs, _ = identify_and_cut_bridges(graph)
    write_subgraphs_to_db(parent_id, subgraphs)

def get_subgraphs_dataframe():
    conn = sqlite3.connect('molecules.db')
    df = pd.read_sql_query('SELECT * FROM subgraphs', conn)
    conn.close()
    return df

def main():
    create_database()

    # Read SMILES strings from a text file
    with open('smiles_fin.txt', 'r') as file:
        initial_smiles = [line.strip() for line in file if line.strip()]

    conn = sqlite3.connect('molecules.db')
    cursor = conn.cursor()
    cursor.executemany('INSERT INTO input_smiles (smiles) VALUES (?)', [(smiles,) for smiles in initial_smiles])
    conn.commit()

    cursor.execute('SELECT id, smiles FROM input_smiles')
    smiles_list = cursor.fetchall()
    conn.close()

    for parent_id, smiles in smiles_list:
        process_smiles(smiles, parent_id)

    # Check subgraph count
    conn = sqlite3.connect('molecules.db')
    cursor = conn.cursor()

    cursor.execute('SELECT COUNT(*) FROM subgraphs')
    subgraphs_count = cursor.fetchone()[0]
    print(f"Number of subgraphs generated: {subgraphs_count}")

    conn.close()

    # Get the subgraphs table as a DataFrame and display it
    df = get_subgraphs_dataframe()
    print(df)

if __name__ == '__main__':
    main()