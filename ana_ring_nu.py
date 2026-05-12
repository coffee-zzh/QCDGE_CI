import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
from tqdm import tqdm

def plot_ringnu_distribution(ring_number_dict,output_figure):

    
    ring_number = {k: len(v) for k, v in ring_number_dict.items()}
    plt.bar(ring_number.keys(), ring_number.values(), color='skyblue')
    for nu in ring_number.keys():
        plt.text(nu, ring_number[nu], str(ring_number[nu]), ha='center', va='bottom',fontsize=18)
    plt.title('Distribution of Ring Numbers',fontsize=18,fontweight='bold')
    plt.xticks(range(max(ring_number_dict.keys()) + 1),fontsize=18)
    plt.yticks(fontsize=18)
    plt.savefig(output_figure)
    plt.close()

def write_to_csv(ring_number_dict,output_csv):

    rows = []
    for nu, mols in ring_number_dict.items():
        for mol in mols:
            rows.append((mol,nu))
    df = pd.DataFrame(rows,columns=['MoleculeID','RingNumber'])
    df.to_csv(output_csv, index=False)
def count_ring_number(csv_file,output_figure,output_csv):

    df = pd.read_csv(csv_file,header=0,index_col='Index')
    mols = df.index.tolist()
    ring_number_dict={0:[],1:[],2:[],3:[],4:[],5:[],6:[],7:[],8:[]}
    ring_number_list = df['RingNumber'].tolist()
    for mol,nu in tqdm(zip(mols,ring_number_list), desc="Counting ring numbers"):
        ring_number_dict[nu].append(mol)

    plot_ringnu_distribution(ring_number_dict,output_figure)
    write_to_csv(ring_number_dict,output_csv)


def main():
    main_path = Path('real_work_dir')
    input_path = main_path / 'final_all.csv'
    ana_ring_path = main_path / 'ana_ring_nu'
    ana_ring_path.mkdir(parents=True, exist_ok=True)
    output_figure = ana_ring_path / 'ring_number.svg'
    output_csv = ana_ring_path / 'ring_number.csv'

    count_ring_number(input_path,output_figure,output_csv)



if __name__ == "__main__":
    main()
