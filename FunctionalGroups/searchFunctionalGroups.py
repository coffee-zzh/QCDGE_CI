import os
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import RDConfig
from rdkit.Chem import FragmentCatalog

def main_searchGroupsImport(smiles):
    # 获取当前脚本文件的绝对路径
    script_path = os.path.abspath(__file__)
    # 获取脚本所在目录的路径
    script_dir = os.path.dirname(script_path)
    # 构建a.dat文件的绝对路径
    fg_file = os.path.join(script_dir, 'FunctionalGroups_for_plot.txt')

    fparams, fcat = getFuncGroupFromFragment(smiles,fg_file)
    groups=[]
    if fcat.GetNumEntries() != 0:
        for i in range(fcat.GetNumEntries()):
            # 向存储器传入分子片段id，获取片段中所包含的官能团编号：fcat.GetEntryFuncGroupIds()
            groups+=list(fcat.GetEntryFuncGroupIds(i))
        groups = set(groups)
        funcnames=[]
        for i in groups:
            funcgroup = fparams.GetFuncGroup(i)
            funcnames.append(f"{funcgroup.GetProp('_Name')}_{i}")
        return groups
    else:
        return None
        # try:
        #     funcnames=process_unsaturated_other_atoms(smiles,fg_file)
        #     return funcnames
        # except:
        #     return None

def fragments(mol,fg_file):
    fparams = FragmentCatalog.FragCatParams(0, 10, fg_file)
    # 传入参数器，创建一个片段存储器,产生的分子片段都会存储在该对象中
    fcat = FragmentCatalog.FragCatalog(fparams)
    # 创建一个片段生成器
    fcgen = FragmentCatalog.FragCatGenerator()
    # 计算分子片段
    fcgen.AddFragsFromMol(mol, fcat)
    #通过存储器查看片段：fcat.GetEntryDescription() 尖括号中的内容：表示与片段相连的官能团
    # print(fcat.GetEntryDescription(11))
    return fparams, fcat

def getFuncGroupFromFragment(smiles,fg_file):
    mol = Chem.MolFromSmiles(smiles)
    num_atoms = mol.GetNumAtoms()
    fparams, fcat = fragments(mol,fg_file)
    # 查看分子片段数量fcat.GetNumEntries()
    if num_atoms <3 or fcat.GetNumEntries()==0:
        index_C = find_unsaturated_carbons(smiles)
        if index_C !=999:
            smiles=smiles[:index_C+1]+'CCCC'+smiles[index_C+1:]
            mol = Chem.MolFromSmiles(smiles)
            fparams, fcat = fragments(mol,fg_file)
            saveImage(mol,'aaaaa')
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
                # 向存储器传入分子片段id，获取片段中所包含的官能团编号：fcat.GetEntryFuncGroupIds()
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


def searchGroupsTest1():
    # smiles='O=CC#N'
    # smiles='O=N(=O)OC1=CN=NO1'
    # smiles='n1cncncnc1'
    # smiles='CC(O)CO'
    # smiles='CCCCC=CC(=O)C'
    #smiles='NC(=O)C(N)=O'
    #smiles='OCCNC1=NC=NO1'
    #smiles = '[NH]C1=C(O)[CH]NC(F)=C1'
    smiles = '[CHCCCC]1C=C(C(=C=O)C=[C]1)O'
    mol = Chem.MolFromSmiles(smiles)
    
    # saveImage(mol,'1.png')
    fg_file = 'FunctionalGroups.txt'
    num_atoms = mol.GetNumAtoms()
    fparams, fcat = getFuncGroupFromFragment(smiles,fg_file)
    groups=[]
    for i in range(fcat.GetNumEntries()):
        # 向存储器传入分子片段id，获取片段中所包含的官能团编号：fcat.GetEntryFuncGroupIds()
        groups+=list(fcat.GetEntryFuncGroupIds(i))
    for i in set(groups):
        # 向参数器传入官能团编号，获取官能团对应的mol对象：fparams.GetFuncGroup()
        funcgroup = fparams.GetFuncGroup(i)
        name = funcgroup.GetProp('_Name')
        print(name+ '    '+Chem.MolToSmarts(funcgroup))
    # 根据id获取片段：fcat.GetEntryDescription()
    # 获取上级片段id：fcat.GetEntryDownIds()

if __name__ == '__main__':
    # searchGroupsTest1()
    main_searchGroupsImport('CCCCC=CC(=O)C')
    # smiles='NCCOC(C)=O'
    # # for i in ['NC(N)=O', 'N#CC#N', 'NC(=N)C#N', 'FC(F)(F)F', 'NC(=O)C(N)=O', 'O=C(C#N)C#N', 'N#CC#CC#N']:
    # a=main_searchGroupsImport(smiles)
    # print(a)


#[C;!$(C=*);!$(CC=O);!a;!$(C[N,O,F])]=N-[N;R0;!$(N=*);!$(NC=O);!$(NC(=O));!a;!$(N[O,F]);!$(NNN)]          
 
