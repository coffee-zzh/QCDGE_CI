import os
import sys
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import RDConfig
from rdkit.Chem import FragmentCatalog


def main_searchGroupsImport(smiles):

    script_path = os.path.abspath(__file__)
    
    script_dir = os.path.dirname(script_path)
    
    fg_file = os.path.join(script_dir, 'FunctionalGroups_for_plot.txt')

    fparams, fcat = getFuncGroupFromFragment(smiles,fg_file)
    groups=[]
    if fcat.GetNumEntries() != 0:
        for i in range(fcat.GetNumEntries()):
            groups+=list(fcat.GetEntryFuncGroupIds(i))
        groups = set(groups)
        funcnames=[]
        for i in groups:
            funcgroup = fparams.GetFuncGroup(i)
            funcnames.append(f"{funcgroup.GetProp('_Name')}_{i}")
        return groups
    else:
        return None


def fragments(mol,fg_file):
    fparams = FragmentCatalog.FragCatParams(0, 10, fg_file)
    fcat = FragmentCatalog.FragCatalog(fparams)
    fcgen = FragmentCatalog.FragCatGenerator()
    try:
        fcgen.AddFragsFromMol(mol, fcat)
    except Exception as e:
        print(f"Warning: AddFragsFromMol Fragmentation failed; treat as an empty result: {e}", file=sys.stderr)
    return fparams, fcat

def getFuncGroupFromFragment(smiles,fg_file):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Warning: MolFromSmiles Parsing failed; treat as an empty result.: {smiles}", file=sys.stderr)
        fparams = FragmentCatalog.FragCatParams(0, 10, fg_file)
        fcat = FragmentCatalog.FragCatalog(fparams)
        return fparams, fcat
    num_atoms = mol.GetNumAtoms()
    fparams, fcat = fragments(mol,fg_file)
    if num_atoms <3 or fcat.GetNumEntries()==0:
        index_C = find_unsaturated_carbons(smiles)
        if index_C !=999:
            smiles=smiles[:index_C+1]+'CCCC'+smiles[index_C+1:]
            mol = Chem.MolFromSmiles(smiles)
            fparams, fcat = fragments(mol,fg_file)

    return fparams, fcat

def find_unsaturated_carbons(smiles):
    mol = Chem.MolFromSmiles(smiles)
    unsaturated_carbons = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetNumImplicitHs()!=0:
            unsaturated_carbons.append(atom.GetIdx())
    if len(unsaturated_carbons)==0:
        last_atom = mol.GetAtomWithIdx(mol.GetNumAtoms() - 1)
        if last_atom.GetNumImplicitHs()!=0:
            unsaturated_carbons.append(atom.GetIdx())

    nu=-1
    Finding=False
    if unsaturated_carbons:
        for c in smiles:
            nu+=1
            if Finding !=True:
                if not c.isalpha():
                    unsaturated_carbons[0]+=1
            if nu == unsaturated_carbons[0]:
                Finding=True
            if Finding == True and nu < len(smiles)-1 and smiles[nu+1] in ['=','(','#']:
                unsaturated_carbons[0]-=1
                break
        return unsaturated_carbons[0]
    else:
        return 999



def find_unsaturated_other_atoms(smiles):
    mol = Chem.MolFromSmiles(smiles)
    _unsaturated = []
    for atom in mol.GetAtoms():
        if atom.GetNumImplicitHs()!=0:
            _unsaturated.append(atom.GetIdx())
    _list=[]
    for i, c in enumerate(smiles):
        if c.isalpha():
            _list.append(i)
    unsaturated_atoms=[]
    for j in _unsaturated:
        unsaturated_atoms.append(_list[j])
    return unsaturated_atoms

def process_unsaturated_other_atoms(smiles,fg_file):
    unsaturated_atoms=find_unsaturated_other_atoms(smiles)
    funcnames_list = []
    for i in unsaturated_atoms:
        smi1=smiles[:i+1]+'(CCCC)'+smiles[i+1:]
        mol = Chem.MolFromSmiles(smi1)
        fparams, fcat = fragments(mol,fg_file)
        groups=[]
        if fcat.GetNumEntries() != 0:
            for j in range(fcat.GetNumEntries()):
                groups+=list(fcat.GetEntryFuncGroupIds(j))
            groups = set(groups)
            funcnames=[]
            for k in groups:
                funcgroup = fparams.GetFuncGroup(k)
                funcnames.append(funcgroup.GetProp('_Name'))
            funcnames_list.append(funcnames)
    funcnames = get_intersection_of_sets(funcnames_list)
    return funcnames

def get_intersection_of_sets(list_of_lists):
    intersection = set(list_of_lists[0])
    for sublist in list_of_lists:
        intersection &= set(sublist)
    return intersection


def saveImage(mol,pngname):
    img=Draw.MolToImage(mol)
    os.makedirs('./images',exist_ok=True)
    img.save(f'./images/{pngname}.png')

def TEST_searchGroupsImport(nu,smiles):
    fg_file = 'FunctionalGroups.txt'
    try:
        mol = Chem.MolFromSmiles(smiles)
        fparams, fcat = getFuncGroupFromFragment(smiles,fg_file)
        groups=[]
        for i in range(fcat.GetNumEntries()):
            groups+=list(fcat.GetEntryFuncGroupIds(i))
        funcnames=[]
        funcgroups=[]
        for i in set(groups):
            funcgroup = fparams.GetFuncGroup(i)
            funcnames.append(funcgroup.GetProp('_Name'))
        return funcnames
    except:
        return None

if __name__ == '__main__':
    main_searchGroupsImport('CCCCC=CC(=O)C')
    
