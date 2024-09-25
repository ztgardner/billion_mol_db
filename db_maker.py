import pandas as pd


def generate_csv(data, n, filename='data.csv'):
    """
    Generates a CSV file with n rows, repeating the provided data.

    Parameters:
        data (dict): A dictionary containing the data for each column. The keys should match the column names.
        n (int): The number of rows to generate.
        filename (str): The name of the output CSV file. Default is 'data.csv'.
    """
    # Create a DataFrame with the provided data
    df = pd.DataFrame([data] * n)

    # Add an 'id' column with unique values
    df['id'] = range(1, n + 1)

    # Save the DataFrame to a CSV file
    df.to_csv(filename, index=False)


# Example usage
ex_data = {
    'd1_smiles': "c5ccc4cc3cc2cc1ccccc1cc2cc3cc4c5",
    'd2_graph': """Vertices: [(0, 'C'), (1, 'C'), (2, 'C'), (3, 'C'), (4, 'C'), (5, 'C'), (6, 'C'), (7, 'C'), (8, 'C'), (9, 'C'),
    (10, 'C'), (11, 'C'), (12, 'C'), (13, 'C'), (14, 'C'), (15, 'C'), (16, 'C'), (17, 'C'), (18, 'C'), (19, 'C'),
    (20, 'C'), (21, 'C'), (22, 'H'), (23, 'H'), (24, 'H'), (25, 'H'), (26, 'H'), (27, 'H'), (28, 'H'), (29, 'H'),
    (30, 'H'), (31, 'H'), (32, 'H'), (33, 'H'), (34, 'H'), (35, 'H')]
    Edges: 0,1: Aromatic, 1,2: Aromatic, 2,3: Aromatic, 3,4: Aromatic, 4,5: Aromatic, 5,6: Aromatic, 6,7: Aromatic,
    7,8: Aromatic, 8,9: Aromatic, 9,10: Aromatic, 10,11: Aromatic, 11,12: Aromatic, 12,13: Aromatic, 13,14: Aromatic,
    14,15: Aromatic, 15,16: Aromatic, 16,17: Aromatic, 17,18: Aromatic, 18,19: Aromatic, 19,20: Aromatic, 20,21: Aromatic,
    14,9: Aromatic, 16,7: Aromatic, 18,5: Aromatic, 20,3: Aromatic, 21,0: Aromatic, 0,22: Single, 1,23: Single, 2,24: Single,
    4,25: Single, 6,26: Single, 8,27: Single, 10,28: Single, 11,29: Single, 12,30: Single, 13,31: Single, 15,32: Single,
    17,33: Single, 19,34: Single, 21,35: Single""",
    'd3_xyz': """C  -6.143599 0.205719  0.247573
          C  -6.052480 -0.968439  -0.493305
          C  -4.800193 -1.480955  -0.835496
          C  -3.622518 -0.824804  -0.440531
          C  -2.354741 -1.327655  -0.776873
          C  -1.177076 -0.670776  -0.381450
          C  0.091088 -1.173738  -0.717870
          C  1.269314 -0.517784  -0.323034
          C  2.537064 -1.021732  -0.660065
          C  3.714925 -0.365939  -0.265322
          C  4.982573 -0.869164  -0.601898
          C  6.143601 -0.205724  -0.202081
          C  6.052479 0.968434  0.538796
          C  4.800193 1.480955  0.880984
          C  3.622517 0.824808  0.486016
          C  2.354741 1.327659  0.822359
          C  1.177076 0.670780  0.426931
          C  -0.091088 1.173741  0.763350
          C  -1.269314 0.517790  0.368514
          C  -2.537065 1.021737  0.705549
          C  -3.714925 0.365942  0.310808
          C  -4.982573 0.869162  0.647388
          H  -7.117682 0.606650  0.515180
          H  -6.955227 -1.486747  -0.805720
          H  -4.749590 -2.400001  -1.415067
          H  -2.283878 -2.246244  -1.356477
          H  0.162380 -2.092402  -1.297526
          H  2.608771 -1.940267  -1.239652
          H  5.074414 -1.785627  -1.180489
          H  7.117685 -0.606655  -0.469686
          H  6.955224 1.486739  0.851219
          H  4.749593 2.399998  1.460559
          H  2.283880 2.246247  1.401967
          H  -0.162380 2.092401  1.343010
          H  -2.608774 1.940269  1.285142
          H  -5.074415 1.785623  1.225984""",
}

n_rows = 3  # Number of rows you want in the CSV
generate_csv(ex_data, n_rows, filename='data_10mil_actual.csv')

ex_data_2 = {
    'd1_smiles': "c5ccc4cc3cc2cc1ccccc1cc2cc3cc4c5",
    'd2_graph': 123,  # Placeholder
    'd3_xyz': 123,  # Placeholder
}

generate_csv(ex_data_2, n_rows, filename='data_10mil_reference.csv')
