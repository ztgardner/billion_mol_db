import networkx as nx
from rdkit import Chem
import matplotlib.pyplot as plt

def smiles_to_graph(smiles):
    """Convert a SMILES string to a graph with atom and bond information."""
    mol = Chem.MolFromSmiles(smiles)
    g = nx.Graph()
    for atom in mol.GetAtoms():
        g.add_node(atom.GetIdx(), element=atom.GetSymbol())
    for bond in mol.GetBonds():
        g.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
    return g

def draw_graph(graph):
    """Draw a graph with atom labels."""
    pos = nx.spring_layout(graph)  # Position nodes using Fruchterman-Reingold force-directed algorithm
    labels = {node: data['element'] for node, data in graph.nodes(data=True)}

    plt.figure(figsize=(10, 7))
    nx.draw(graph, pos, labels=labels, with_labels=True, node_color='lightblue', node_size=3000, font_size=16)
    plt.show()

# Example usage
smiles_string = "Oc1ccccc1c2cccs2"
graph = smiles_to_graph(smiles_string)
draw_graph(graph)