import random

def select_random_smiles(input_file, output_file, n):
    # Read SMILES strings from the input file
    with open(input_file, 'r') as file:
        smiles_list = [line.strip() for line in file.readlines()]

    # Check if n is larger than the number of available SMILES
    if n > len(smiles_list):
        raise ValueError("n is larger than the number of lines in the input file")

    # Randomly select n SMILES strings
    selected_smiles = random.sample(smiles_list, n)

    # Write the selected SMILES strings to the output file
    with open(output_file, 'w') as file:
        for smiles in selected_smiles:
            file.write(smiles + '\n')

if __name__ == "__main__":
    input_file = 'smiles.txt'
    n = 10 # Number of lines to randomly select
    output_file = f'smiles_{n}.txt'

    select_random_smiles(input_file, output_file, n)