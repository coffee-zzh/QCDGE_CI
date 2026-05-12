
import numpy as np
from pathlib import Path
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from ase.data import chemical_symbols
from scipy import stats
import h5py
import tempfile
import os
import json
from tqdm import tqdm
import time

try:
    import openbabel as ob
except ImportError:
    ob = None


def write_xyz(atoms, coords, filename):
    
    with open(filename, 'w') as f:
        f.write(f"{len(atoms)}\n")
        f.write(f"Ground State\n")
        for i, (atom, coord) in enumerate(zip(atoms, coords)):
            f.write(f"{atom} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n")


def get_bond_pairs_from_xyz(ground_xyz):
    
    if ob is None:
        print("Warning: OpenBabel not found, using RDKit fallback (only single bonds)")
        mol = Chem.MolFromXYZFile(ground_xyz)
        if mol is None:
            return []
        rdDetermineBonds.DetermineConnectivity(mol)
        bonds = []
        for bond in mol.GetBonds():
            idx1 = bond.GetBeginAtomIdx()
            idx2 = bond.GetEndAtomIdx()
            elem1 = bond.GetBeginAtom().GetSymbol()
            elem2 = bond.GetEndAtom().GetSymbol()
            if idx1 > idx2:
                idx1, idx2 = idx2, idx1
                elem1, elem2 = elem2, elem1
            bond_type = f"{elem1}-{elem2}" if elem1 < elem2 else f"{elem2}-{elem1}"
            bonds.append((idx1, idx2, elem1, elem2, bond_type, 1.0))
        return sorted(bonds, key=lambda x: (x[0], x[1]))

   
    obconversion = ob.OBConversion()
    obconversion.SetInAndOutFormats('xyz', 'mol')

    obmol = ob.OBMol()
    obconversion.ReadFile(obmol, str(ground_xyz))
    obmol.ConnectTheDots()
    obmol.PerceiveBondOrders()
    obmol.SetAromaticPerceived()

    bonds = []
    for bond in ob.OBMolBondIter(obmol):
        idx1 = bond.GetBeginAtomIdx() - 1  
        idx2 = bond.GetEndAtomIdx() - 1
        elem1 = ob.GetSymbol(bond.GetBeginAtom().GetAtomicNum())
        elem2 = ob.GetSymbol(bond.GetEndAtom().GetAtomicNum())

        
        if bond.IsAromatic():
            bond_order = 1.5 
        else:
            bond_order = float(bond.GetBondOrder())  

        if idx1 > idx2:
            idx1, idx2 = idx2, idx1
            elem1, elem2 = elem2, elem1

        bond_type = f"{elem1}-{elem2}" if elem1 < elem2 else f"{elem2}-{elem1}"
        bonds.append((idx1, idx2, elem1, elem2, bond_type, bond_order))

    return sorted(bonds, key=lambda x: (x[0], x[1]))


def calc_bond_lengths(coords, bond_pairs):
    
    lengths = {}
    for idx1, idx2, elem1, elem2, bond_type, bond_order in bond_pairs:
        dist = np.linalg.norm(coords[idx1] - coords[idx2])
        lengths[(idx1, idx2)] = dist
    return lengths


def plot_bond_type_comparison(bond_data, outpath, bw_method=None):

    for bond_type in sorted(bond_data.keys()):
        all_ground = []
        all_ci = []
        for bond_order in bond_data[bond_type].keys():
            all_ground.extend(bond_data[bond_type][bond_order]['ground'])
            all_ci.extend(bond_data[bond_type][bond_order]['ci'])

        if len(all_ground) < 1:
            continue

        has_h = 'H' in bond_type
        has_f = 'F' in bond_type
        if has_h:
            x_min, x_max = 0.7, 1.5 
            current_method = 0.5
        elif has_f:
            x_min, x_max = 1, 2
            current_method = 0.5
        else:
            x_min, x_max = 1, 2
            current_method = bw_method

        x_full = np.linspace(x_min, x_max, 500)
        kde_ground = stats.gaussian_kde(all_ground, bw_method=current_method)
        kde_ci = stats.gaussian_kde(all_ci, bw_method=current_method)
        y_ground_full = kde_ground(x_full)
        y_ci_full = kde_ci(x_full)
        global_max = max(y_ground_full.max(), y_ci_full.max())

        
        mask_short = (x_full >= x_min) & (x_full <= x_max)
        x_plot_short = x_full[mask_short]
        y_ground_plot_short = y_ground_full[mask_short] / global_max
        y_ci_plot_short = y_ci_full[mask_short] / global_max

        fig, ax1 = plt.subplots(figsize=(18, 10 * 18 / 16))
        ax1.fill_between(x_plot_short, y_ground_plot_short, alpha=0.6, color="#2ca8c5", label='GS')
        ax1.fill_between(x_plot_short, y_ci_plot_short, alpha=0.6, color="#98d6b9", label='CI')
        ax1.plot(x_plot_short, y_ground_plot_short, color="#2ca8c5", linewidth=2.5)
        ax1.plot(x_plot_short, y_ci_plot_short, color="#98d6b9", linewidth=2.5)

        ax1.set_xlabel('Bond Length (Å)', fontsize=50)
        ax1.set_ylabel('Probability density', fontsize=50)
        ax1.set_title(f'{bond_type} Distribution (Short Bonds)', fontsize=50)
        ax1.set_ylim(bottom=0)
        ax1.set_xlim(x_min, x_max)
        ax1.tick_params(labelsize=30)
        ax1.legend(fontsize=50)
        plt.tight_layout()
        plt.savefig(outpath / f"{bond_type.replace('-','_')}_kde_comparison_short.svg", dpi=500)
        plt.close()


def plot_bond_order_comparison(bond_data, outpath, bw_method=None):
   
    for bond_type in sorted(bond_data.keys()):
        for bond_order in sorted(bond_data[bond_type].keys()):
            data_ground = bond_data[bond_type][bond_order]['ground']
            data_ci = bond_data[bond_type][bond_order]['ci']

            if len(data_ground) <= 1:
                continue

            has_h = 'H' in bond_type
            has_f = 'F' in bond_type
            if has_h or has_f:
                continue
            x_min, x_max = 1, 2
            if bond_order == 1.5:
                current_method = 0.5
            else:
                current_method = bw_method

            x_full = np.linspace(x_min, x_max, 500)
            kde_ground = stats.gaussian_kde(data_ground, bw_method=current_method)
            kde_ci = stats.gaussian_kde(data_ci, bw_method=current_method)
            y_ground_full = kde_ground(x_full)
            y_ci_full = kde_ci(x_full)
            global_max = max(y_ground_full.max(),y_ci_full.max())

            fig, ax1 = plt.subplots(figsize=(18, 10 * 18 / 16))
            mask = (x_full >= x_min) & (x_full <= x_max)
            x_plot = x_full[mask]
            y_ground_plot = y_ground_full[mask] / global_max
            y_ci_plot = y_ci_full[mask] / global_max

            ax1.fill_between(x_plot, y_ground_plot, alpha=0.6, color="#2ca8c5", label='GS')
            ax1.fill_between(x_plot, y_ci_plot, alpha=0.6, color="#98d6b9", label='CI')
            ax1.plot(x_plot, y_ground_plot, color="#2ca8c5", linewidth=2.5)
            ax1.plot(x_plot, y_ci_plot, color="#98d6b9", linewidth=2.5)

            bo_label = "aromatic" if bond_order == 1.5 else str(bond_order)
            ax1.set_xlabel('Bond Length (Å)', fontsize=50)
            ax1.set_ylabel('Probability density', fontsize=50)
            ax1.set_title(f'{bond_type} (bond order {bo_label}) Distribution', fontsize=50)
            ax1.set_xlim(x_min, x_max)
            ax1.tick_params(labelsize=30)
            ax1.legend(fontsize=50)
            plt.tight_layout()
            bo_str = "aromatic" if bond_order == 1.5 else str(bond_order)
            plt.savefig(outpath / f"{bond_type.replace('-','_')}_order_{bo_str}_kde_comparison_short.svg", dpi=500)
            plt.close()
            


def analyze_pair(hdf5_path):

    bond_data = defaultdict(lambda: defaultdict(lambda: {'ground': [], 'ci': []}))
    all_data = []

    with h5py.File(hdf5_path, 'r') as f:
        groups = list(f.keys())
        for group in tqdm(groups, desc="Analyzing molecules", unit="molecule"):
            print(group)
            group_data = f[group]
            labels = group_data.attrs['labels'][()]
            atoms = [chemical_symbols[i] for i in labels]
            ground_coord = group_data['Ground']["coordinates"][()]
            ci_coord = group_data['CI']["coordinates"][()]

            with tempfile.NamedTemporaryFile(delete=False, suffix=".xyz") as temp_file:
                xyz_file = temp_file.name
                write_xyz(atoms, ground_coord, xyz_file)
            bond_pairs = get_bond_pairs_from_xyz(xyz_file)
            os.remove(xyz_file)

            ground_lengths = calc_bond_lengths(ground_coord, bond_pairs)
            ci_lengths = calc_bond_lengths(ci_coord, bond_pairs)

            

            for idx1, idx2, elem1, elem2, bond_type, bond_order in bond_pairs:
                gl = ground_lengths[(idx1, idx2)]
                cl = ci_lengths[(idx1, idx2)]
                
                    
                diff = cl - gl
                bond_data[bond_type][bond_order]['ground'].append(gl)
                bond_data[bond_type][bond_order]['ci'].append(cl)
                all_data.append((idx1, idx2, elem1, elem2, bond_type, bond_order, gl, cl, diff))

    print("="*100)
    print(f"{'Bond Type':<10} {'Order':<10} {'Count':>8} {'Mean(Ground)':>12} {'Mean(CI)':>12} {'Diff':>12} {'Ground_bond_range':>25} {'CI_bond_range':>25}")
    print("="*100)

    for bt in sorted(bond_data.keys()):
        for bo in sorted(bond_data[bt].keys()):
            g_list = bond_data[bt][bo]['ground']
            c_list = bond_data[bt][bo]['ci']
            bo_label = "aromatic" if bo == 1.5 else str(bo)
            if len(g_list) > 0:
                print(f"{bt:<10} {bo_label:<10} {len(g_list):>8} "
                      f"{np.mean(g_list):>12.4f} {np.mean(c_list):>12.4f} {np.mean(c_list)-np.mean(g_list):>12.4f} "
                      f"{min(g_list):>12.4f} - {max(g_list):>12.4f} {min(c_list):>12.4f} - {max(c_list):>12.4f}")

    print("="*100)

    return bond_data, all_data


def main():
    main_path = Path("real_work_dir")
    ana_path = main_path / "ana_bond_types"
    start_time = time.time()
    ana_path.mkdir(parents=True, exist_ok=True)
    hdf5_path = main_path / "Final_property.hdf5"
    bond_data, all_data = analyze_pair(hdf5_path)
    bw_method = 0.25

    plot_bond_type_comparison(bond_data, ana_path, bw_method=bw_method)
    plot_bond_order_comparison(bond_data, ana_path, bw_method=bw_method)
    with open(ana_path / "all_bond_data.json", 'w') as f:
        json.dump(all_data, f, indent=4)
    with open(ana_path / "bond_data.json",'w') as f:
        json.dump(bond_data, f, indent=4)
    
    end_time = time.time()
    print(f"Total time: {end_time - start_time:.2f} seconds")

if __name__ == '__main__':
    main()

