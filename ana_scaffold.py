import pandas as pd
import pickle
from collections import Counter
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Scaffolds import MurckoScaffold
from pathlib import Path


def analyzeScaffold(smiles_list, inchi_list,ana_scaffold_path):
    mols=[]
    for smiles, inchi in zip(smiles_list, inchi_list):

        mol = None
        if smiles != '1':
            mol = Chem.MolFromSmiles(smiles)
        if mol is None:  
            mol = Chem.MolFromInchi(inchi)  

        mols.append(mol)

    mol_scaffolds = [MurckoScaffold.GetScaffoldForMol(mol) for mol in mols if mol]
    
    if label == 'C':
        grafh_scaffolds = [MurckoScaffold.MakeScaffoldGeneric(s) for s in mol_scaffolds]
    
    else:
        grafh_scaffolds = mol_scaffolds
    scaffold_smiles = [Chem.MolToSmiles(scaffold) for scaffold in grafh_scaffolds if scaffold != None]
 
    lists=[]
    for nu,(i,j) in enumerate(zip(mols, scaffold_smiles)):
        if i. GetRingInfo().NumRings() == 0:
            if j == '':
                lists.append(smiles_list.iloc[nu])
    with open(ana_scaffold_path / f'none_scaffoldsmiles_{label}.log', 'w') as file:
        for item in lists:
            file.write(f"{item}\n")
    counter=Counter(scaffold_smiles)
    with open(ana_scaffold_path / f'MurckoScaffold_all_{label}.pkl', 'wb') as file:
        pickle.dump(counter, file)

    sorted_dict = dict(sorted(counter.items(), key=lambda item: item[1], reverse=True)[:21])
    for key, value in sorted_dict.items():
        print(f'{key}: {value}')
    len_functional_groups_count = len(counter)
    print(f"Length of the loaded object: {len_functional_groups_count}")


def analyzeScaffold_main(input_file,ana_scaffold_path):
    smiles_list = pd.read_csv(input_file, index_col=0)[col_label1]
    inchi_list = pd.read_csv(input_file, index_col=0)[col_label2]
    analyzeScaffold(smiles_list, inchi_list,ana_scaffold_path)

def get_scaffold_fig(df,ana_scaffold_path):

    sorted_dict = sorted(df.items(), key=lambda x: x[1], reverse=True)[:21]

    for key, value in sorted_dict:
        print(f"{key}: {value}")

    smis = [i[0] for i in sorted_dict]
    print(smis)
    molecules = [Chem.MolFromSmiles(smiles) for smiles in smis]

    # Generate molecule image
    img = Draw.MolsToGridImage(molecules, molsPerRow=4, subImgSize=(200,200))
    img.save(ana_scaffold_path / f'scaffolds_img_{label}.png')


def pie_plot(data,ana_scaffold_path):
    
    #Filtered data which value > 2000, and calculate the sum of the rest
    filtered_data = {k: v for k, v in data.items() if v > 2000}
    filtered_data["Other"] = sum(v for k, v in data.items() if v <= 2000)

    # Prepare data
    labels = filtered_data.keys()
    sizes = filtered_data.values()

    # draw pie chart
    plt.figure(figsize=(10, 10))
    plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=140, textprops={'fontsize': 8}, pctdistance=0.85)
    plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

    plt.title('Pie Chart of Categories with Values > 2000 and Others')
    plt.tight_layout()
    plt.savefig(ana_scaffold_path / f'scaffold_pie_{label}.png')

def plot_main(ana_scaffold_path):
    with open(ana_scaffold_path / f'MurckoScaffold_all_{label}.pkl', 'rb') as file:
        df = pickle.load(file)
    pie_plot(df,ana_scaffold_path)
    get_scaffold_fig(df,ana_scaffold_path)

if __name__ == '__main__':
    main_path = Path('real_work_dir')
    ana_scaffold_path = main_path / "ana_scaffold"
    ana_scaffold_path.mkdir(parents=True, exist_ok=True)
    global label
    choose = int(input('(1) Normal or (2) treat all atoms as C and all bonds as single type\n'))

    label='atom' if choose == 1 else 'C'
    
    smi_file = main_path / 'final_all.csv'

    global col_label1, col_label2
    col_label1 = 'Smiles_rdkit'
    col_label2 = 'Inchi_rdkit'

    analyzeScaffold_main(smi_file,ana_scaffold_path)
    plot_main(ana_scaffold_path)

