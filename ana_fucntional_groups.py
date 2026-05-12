import pickle
import pandas as pd
import subprocess as sub
from collections import Counter
from pathlib import Path

from FunctionalGroups.searchFunctionalGroups import main_searchGroupsImport as search_FG
from rdkit.Chem import FragmentCatalog

from rdkit import Chem
from rdkit.Chem import Draw

def analyzeFGs_main(input_file):
    counter = Counter()
    none_list = []
    df = pd.read_csv(input_file, index_col=0)['Smiles_rdkit']
    error_list= []

    for smi in df:
        try:
            fgs_id = search_FG(smi)
        except:
            error_list.append(smi)
            print(f'Error: {len(error_list)}')
            continue
        if fgs_id == None:
            none_list.append(smi)
        counter.update(fgs_id)
    all_dict = {
    "count_id":counter,
    "list": none_list
    }
    with open(f'{ana_fungroup}/FunctionalGroups_count_rdkit.pkl', 'wb') as file:
        pickle.dump(all_dict, file)


def draw_img_of_top(ana_fungroup):
    with open(f'{ana_fungroup}/FunctionalGroups_count_rdkit.pkl', 'rb') as file:
        df = pickle.load(file)
    fg_file=f'{main_path}/FunctionalGroups_for_plot.txt'
    fparams = FragmentCatalog.FragCatParams(1,6,fg_file)
    mols=[]
    top=25
    sorted_dict = sorted(df['count_id'].items(), key=lambda x: x[1], reverse=True)[:top]
    for key, value in sorted_dict:
        funcgroup = fparams.GetFuncGroup(key)
        name=f"{funcgroup.GetProp('_Name')}_{key}"
        print(f"{name}: {value}")
        mols.append(funcgroup)
    img=Draw.MolsToGridImage(mols,molsPerRow=6)
    img.save(f'{ana_fungroup}/fgs_of_top_{top}.png')
    print(sorted(df['count_id'].keys()))
    print(len(df['count_id']))
    print(df['list'])


def checkmol_main(input_file):
    df = pd.read_csv(input_file, index_col=0)['Smiles_rdkit']
    smiles_to_mol(df)
    run_checkmol(df)
    

def smiles_to_mol(df):
    for nu,smi in enumerate(df):
        mol = Chem.MolFromSmiles(smi)
        mol_block = Chem.MolToMolBlock(mol)
        mol_file_path = f'./FunctionalGroups/checkmol/mols/index_{nu}.mol'
        with open(mol_file_path, 'w') as f:
            f.write(mol_block)

def run_checkmol(df):
    counter = Counter()
    for nu,_ in enumerate(df):
        result = sub.run(f'./FunctionalGroups/checkmol/checkmol-0.5b-linux-x86_64 -e ./FunctionalGroups/checkmol/mols/index_{nu}.mol',shell=True, capture_output=True, text=True)
        out = result.stdout.strip().split('\n')
        counter.update(out)
    with open(f'{ana_fungroup}/FunctionalGroups_count_checkmol.pkl', 'wb') as file:
        pickle.dump(counter, file)

def checkmol_top(ana_fungroup):
    with open(f'{ana_fungroup}/FunctionalGroups_count_checkmol.pkl', 'rb') as file:
        df = pickle.load(file)
    print("checkmol result:")
    print(len(df))

if __name__ == '__main__':
    main_path = Path("real_work_dir")
    ana_fungroup = main_path / 'ana_functional_groups'
    ana_fungroup.mkdir(parents=True, exist_ok=True)
    smi_file = main_path / 'final_all.csv'
    import time
    start_time = time.time()
    analyzeFGs_main(smi_file)
    # draw_img_of_top(ana_fungroup)
    # checkmol_main(smi_file)
    # checkmol_top(ana_fungroup)
    end_time = time.time()
    print(f"Total time: {end_time - start_time} seconds")
    
