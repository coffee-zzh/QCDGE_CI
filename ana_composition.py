import os
import re
import numpy as np
import pandas as pd
import json
from rdkit import Chem
from collections import defaultdict
from sub_ana_classify_mol import classify_main
from pathlib import Path


def test(df,col_label1):
    for i in df[col_label1]:
        print(i)
        a= Chem.MolFromSmiles(i)



def count_heavy_atoms(df,col_label1, col_label2):
    def calculate_heavy_atom_count(row, smiles_col, inchi_col):
        if row[smiles_col] != '1':
            mol = Chem.MolFromSmiles(row[smiles_col])
        else:
            mol = None

        if mol is None:
            mol = Chem.MolFromInchi(row[inchi_col])

        if mol is not None:
            return mol.GetNumHeavyAtoms()
        else:
            return None

    df['HeavyAtomCount'] = df.apply(calculate_heavy_atom_count, axis=1, smiles_col=col_label1, inchi_col=col_label2)
    return df

def count_cof_molecules(count_dict, smiles_list, inchi_list):
    for smiles, inchi in zip(smiles_list, inchi_list):
        finding = False
        if smiles != '1':
            mol = Chem.MolFromSmiles(smiles)
        else:
            mol = None

        if mol is None:
            mol = Chem.MolFromInchi(inchi)
        if mol is not None:
            atom_symbols = set(atom.GetSymbol() for atom in mol.GetAtoms())
            atom_symbols.discard('H')
            for ele, count in count_dict.items():
                elements = re.findall(r'[A-Z][a-z]*', ele)
                if atom_symbols == set(elements):
                    finding=True
                    count_dict[ele] += 1
                    break
        else:
            print('None!!!!')
        if finding == False:
            print(atom_symbols)
            print(f'SMILES:{smiles}')
    return count_dict


def count_rings_for_each_mol(df, smiles_col, inchi_col):
    def count_rings(row):
        if row[smiles_col] != '1':
            mol = Chem.MolFromSmiles(row[smiles_col])
        else:
            mol = None
        if mol is None:
            mol = Chem.MolFromInchi(row[inchi_col])

        if mol is not None:
            return mol.GetRingInfo().NumRings()
        else:
            return None
    df['RingNumber'] = df.apply(count_rings, axis=1)
    return df


def process_seperately_according_to_HeavyAtomCount(ana_composition_path,df,input_file,col_label1,col_label2,element_str):
    output_prefix, output_extension = os.path.splitext(os.path.basename(input_file))

    split_data = df.groupby('HeavyAtomCount')
    element_compostion = {}

    for group, data in split_data:
        count_dict = dict.fromkeys(element_str, 0)
        element_compostion[f'atom_{group}']=count_cof_molecules(count_dict, data[col_label1].tolist(), data[col_label2].tolist())
        data.reset_index(drop=True, inplace=True)
        data.index+=1
        output_file = ana_composition_path / f"{output_prefix}_{group}_atoms{output_extension}"
        data.to_csv(output_file, index=True)
    save_element_compostion(ana_composition_path, element_compostion)

def save_element_compostion(ana_composition_path, element_compostion):
    with open(ana_composition_path / "1_element_compostion.json", 'w') as json_file:
        json.dump(element_compostion, json_file)

def count_main(main_path,input_file, output_file,element_str):
    col_label1 = 'Smiles_rdkit'
    col_label2 = 'Inchi_rdkit'
    df = pd.read_csv(input_file, low_memory=False)

    df = count_heavy_atoms(df,col_label1, col_label2)
    df = count_rings_for_each_mol(df,col_label1, col_label2)
    df = classify_main(df, col_label1, col_label2)
    ana_composition_path = main_path / "ana_composition"
    process_seperately_according_to_HeavyAtomCount(ana_composition_path,df,input_file,col_label1, col_label2,element_str)
    df.to_csv(output_file)


if __name__ == '__main__':
    main_path = Path("real_work_dir")
    ana_composition_path = main_path / "ana_composition"
    ana_composition_path.mkdir(parents=True, exist_ok=True)
    smi_file = main_path / "Final_Smiles_Inchi.csv"
    output_file = main_path / "final_all.csv"
    element_str=['C','N','O','F','NO','OF','NF','CN','CO','CF', 'NOF', 'CNO','CNF','COF','CNOF']

    count_main(main_path,smi_file,output_file,element_str)
