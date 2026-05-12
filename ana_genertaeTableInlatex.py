import pickle
from rdkit import Chem
from rdkit.Chem import Draw
import os 
from pathlib import Path

def generateTableInLatex(main_path, file_path,label):
    with open(file_path, 'rb') as file:

        functional_groups_count = pickle.load(file)

    total_nu=20

    sorted_dict = dict(sorted(functional_groups_count.items(), key=lambda item: item[1], reverse=True)[:21])

    keys=[]
    values=[]
    for key, value in sorted_dict.items():
        if key != '':
            keys.append(key)
            values.append(value)
        if len(keys) == total_nu:
            break

    molecules = [Chem.MolFromSmiles(smiles) for smiles in keys]

    os.makedirs(main_path / f"scaffold_{label}", exist_ok=True)

    for i, mol in enumerate(molecules):
        print(i+1)
        img = Draw.MolToImage(mol, size=(200, 200))
        img.save(main_path / f"scaffold_{label}/{label}{i+1}.eps")


    a=''
    for i in range(int(total_nu/2)):
        b=f'{i+1} &\\includegraphics[width=2.2cm,height=1.8cm,keepaspectratio]{{{label}{i+1}.eps}} & {keys[i]} & {values[i]} & {int(total_nu/2+1+i)} &\\includegraphics[width=2cm,height=1.8cm,keepaspectratio]{{{label}{int(total_nu/2+1+i)}.eps}} & {keys[int(total_nu/2+i)]} & {values[int(total_nu/2+i)]} \\\\ \\hline\n'
        a+=b
    print(a)

if __name__ == '__main__':
    main_path = Path('real_work_dir')
    ana_scaffold_path = main_path / "ana_scaffold"
    file_path = ana_scaffold_path / 'MurckoScaffold_all_C.pkl'
    label='C'
    generateTableInLatex(ana_scaffold_path, file_path,label)
