import rdkit
from rdkit import Chem
import pandas as pd
import networkx as nx

def smiles_2_rings(smiles):
    """
    Converts a SMILES string to a DataFrame with atom indices and backbone identification.
    """
    # Convert SMILES to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string.")

    # Create a graph from molecule
    G = nx.Graph()
    for bond in mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    # Identify ring atoms using NetworkX cycle detection
    ring_atoms = set()
    cycles = nx.cycle_basis(G)
    for cycle in cycles:
        ring_atoms.update(cycle)

    # Create DataFrame
    data = []
    for atom in mol.GetAtoms():
        atom_idx = atom.GetIdx()
        backbone = "yes" if atom_idx in ring_atoms else "no"
        data.append({"atom": atom_idx, "in_ring": backbone})

    df = pd.DataFrame(data)
    return df

# Example Usage
smiles = "CC(C)[Si](C#CC1=C2C=C3C=CC=CC3=CC2=C(C4=CC5=CC=CC=C5C=C41)C#C[Si](C(C)C)(C(C)C)C(C)C)(C(C)C)C(C)C"  # Example cyclobutane with a sidechain
df = smiles_2_rings(smiles)
print(df)