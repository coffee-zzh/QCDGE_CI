import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm

class CompoundTypeanalyze:
    
    Compound_TYPES = ['carboacyclic', 'heteroacyclic', 'carbocycles', 'heterocycles', 'heteroaromatics', 'fused carbocycles', 'fused heterocycles', 'aromatics']
    
    def __init__(self, compound_type_file):
        self.compound_type_file = compound_type_file  
        self.df, self.mols = self._load_csv(compound_type_file)
    
    def _load_csv(self, compound_type_file):
        df = pd.read_csv(compound_type_file,index_col=0)
        mols = list(df.index)
        return df,mols
    
    def get_mol_type(self, file_path):
        compound_type = {ctype: [] for ctype in self.Compound_TYPES}
        # lines = json.loads(file_path.read_text())
        mols = self.mols
        for mol_id in tqdm(mols, desc='Processing files'):
            ctype = self.df.loc[mol_id, 'CompoundType']
            if ctype in compound_type:
                compound_type[ctype].append(mol_id)
            else:
                compound_type[ctype] = [mol_id]
        return compound_type
    def write_to_csv(self, compound_type, output_file):
        
        rows = []
        for ctype, mol_ids in compound_type.items():
            for mol_id in mol_ids:
                rows.append((mol_id, ctype))
        df = pd.DataFrame(rows, columns=['MoleculeID','CompoundType'])
        df.to_csv(output_file,index=False)

    def draw_distribution_picture(self, file_path, output_figure, output_csv):

        compound_type = self.get_mol_type(file_path)
        self.write_to_csv(compound_type, output_csv)
        ctype_count = {ctype: len(mol_ids) for ctype, mol_ids in compound_type.items()}
        sorted_items = sorted(ctype_count.items(), key=lambda x: x[1], reverse=True)
        sorted_categories = [item[0] for item in sorted_items]
        sorted_counts = [item[1] for item in sorted_items]
        

        plt.figure(figsize=(18, 10*18/16))
        bars = plt.barh(sorted_categories, sorted_counts, color='#96B6D8')
        for bar, count in zip(bars, sorted_counts):
            plt.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height()/2.,
                    f'{count}', ha='left', va='center', fontsize=9)
        
        
        plt.xlabel('Count of molecules')
        plt.ylabel('Molecular categories')
        plt.title(f'Distribution of Molecular categories')
        plt.savefig(output_figure)
        plt.close()
    


def main():
    main_path = Path('real_work_dir')
    ana_compound_type_path = main_path / "ana_compound_type"
    ana_compound_type_path.mkdir(parents=True, exist_ok=True)
    input_path = main_path / 'final_all.csv'
    output_figure = ana_compound_type_path / 'CompoundTypeDistribution.svg'
    output_csv = ana_compound_type_path / 'CompoundTypeDistribution.csv'
    
    
    analyzer = CompoundTypeanalyze(input_path)
    analyzer.draw_distribution_picture(input_path, output_figure, output_csv)
    


if __name__ == "__main__":
    main()

    
