#!/usr/bin/env python3
"""
Molecular Properties Calculator for CSV Files

This script processes a CSV file containing SMILES sequences, calculates various 
molecular properties for each sequence, and saves the results to a new CSV file.

Properties calculated:
- Molecular Weight (MW)
- Hydrophobicity (LogP)
- Number of Hydrogen Bond Acceptors (HBA)
- Number of Hydrogen Bond Donors (HBD)
- Number of Rotatable Bonds (RB)
- Number of Aromatic Rings (arR)

Usage:
    python create_test_data.py input.csv COLUMNNAME output.csv
    
Example:
    python create_test_data.py testdata/input.csv SMILES original-java/data.csv

"""

import sys
import csv
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
from rdkit import RDLogger

# Suppress RDKit warning messages
RDLogger.DisableLog('rdApp.*')


def calculate_properties(smiles):
    """
    Calculate molecular properties for a given SMILES string.
    
    Args:
        smiles (str): The SMILES representation of the molecule
        
    Returns:
        dict: Dictionary containing the calculated properties
    """
    # Convert SMILES to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        return {
            "MW": None,
            "LogP": None,
            "HBA": None,
            "HBD": None,
            "RB": None,
            "arR": None
        }
    
    # Calculate properties
    properties = {
        "MW": round(Descriptors.MolWt(mol), 2),
        "LogP": round(Descriptors.MolLogP(mol), 2),
        "HBA": Lipinski.NumHAcceptors(mol),
        "HBD": Lipinski.NumHDonors(mol),
        "RB": Descriptors.NumRotatableBonds(mol),
        "arR": Chem.Lipinski.NumAromaticRings(mol)
    }
    
    return properties


def process_csv(input_file, smiles_column, output_file):
    """
    Process a CSV file, calculate properties for each SMILES, and save to a new file.
    
    Args:
        input_file (str): Path to input CSV file
        smiles_column (str): Name of the column containing SMILES strings
        output_file (str): Path to output CSV file
    """
    try:
        # Read the CSV file using pandas
        df = pd.read_csv(input_file)
        
        # Check if the SMILES column exists
        if smiles_column not in df.columns:
            print(f"Error: Column '{smiles_column}' not found in the input file")
            return
        
        # Initialize new columns with None values
        df["MW"] = None
        df["LogP"] = None
        df["HBA"] = None
        df["HBD"] = None
        df["RB"] = None
        df["arR"] = None
        
        # Counter for progress reporting
        total_compounds = len(df)
        print(f"Processing {total_compounds} compounds...")
        
        # Process each SMILES
        for idx, row in df.iterrows():
            if idx % 100 == 0 and idx > 0:
                print(f"Processed {idx}/{total_compounds} compounds")
                
            smiles = row[smiles_column]
            
            # Skip empty SMILES
            if pd.isna(smiles) or smiles.strip() == "":
                continue
                
            # Calculate properties
            try:
                properties = calculate_properties(smiles)
                
                # Update the dataframe
                for prop, value in properties.items():
                    df.at[idx, prop] = value
                    
            except Exception as e:
                print(f"Error processing SMILES '{smiles}': {str(e)}")
        
        # Save the results to a new CSV file
        df.to_csv(output_file, index=False, sep ='\t')
        print(f"Results saved to {output_file}")
        output_file_txt = output_file.replace('.csv', '.txt')
        df.to_csv(output_file_txt, index=None, sep='\t', mode='a')
        
    except Exception as e:
        print(f"Error processing CSV file: {str(e)}")


def main():
    """Main function to handle command line arguments."""
    if len(sys.argv) != 4:
        print("Usage: python mol_properties_csv.py input.csv smiles_column_name output.csv")
        sys.exit(1)
    
    input_file = sys.argv[1]
    smiles_column = sys.argv[2]
    output_file = sys.argv[3]
    
    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found")
        sys.exit(1)
    
    # Check if output file already exists
    if os.path.exists(output_file):
        response = input(f"Output file '{output_file}' already exists. Overwrite? (y/n): ")
        if response.lower() != 'y':
            print("Operation cancelled")
            sys.exit(0)
    
    process_csv(input_file, smiles_column, output_file)


if __name__ == "__main__":
    main()