import h5py
from pathlib import Path
from tqdm import tqdm
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy import stats

def classify_compound_type(category_df,mol_list):
    Compound_TYPES = {'carboacyclic':[], 'heteroacyclic':[], 'carbocycles':[], 'heterocycles':[], 'heteroaromatics':[], 'fused carbocycles':[], 'fused heterocycles':[], 'aromatics':[]}

    for mol in mol_list:
        compound_type = category_df.loc[mol,'CompoundType']
        Compound_TYPES[compound_type].append(mol)
    return Compound_TYPES

def read_hdf5_file(hdf5_path, energy_file):
    mol_list = []
    CI_energies = []
    VEE_energies = []
    with h5py.File(hdf5_path, 'r') as f:
        groups = list(f.keys())
        for group in tqdm(groups, desc="Processing groups"):
            group_data = f[group]
            e_CI = group_data['CI']['energy'][()]
            e_Ground = group_data['Ground']['energy'][()]
            CI_energy = e_CI - e_Ground
            CI_energies.append(CI_energy)
            mol_list.append(group)
            VEE = group_data['Ground']['VEE'][()]
            VEE_energies.append(VEE)
    delta_energy_df = pd.DataFrame({'Molecule':mol_list,'CI_energy':CI_energies,'VEE':VEE_energies})
    delta_energy_df.to_csv(energy_file, index=False)

    return delta_energy_df

def plot_all_energy(delta_energy_df,out_path):
    """Draw energy distribution scatter plot for all molecules"""
    fig = plt.figure(figsize=(18, 10*18/16))
    ax1 = fig.add_subplot(111)
    x_vals = np.linspace(0, 12, 400)
    VEE_energies = delta_energy_df['VEE']
    CI_energies = delta_energy_df['CI_energy']
    kde_vee = stats.gaussian_kde(VEE_energies)
    kde_ci = stats.gaussian_kde(CI_energies)
    y_vee = kde_vee(x_vals)
    y_ci = kde_ci(x_vals)
    # Normalize to [0, 1]
    y_max = max(y_vee.max(), y_ci.max())
    y_vee = y_vee / y_max
    y_ci = y_ci / y_max
    ax1.fill_between(x_vals, y_vee, alpha=0.6, color="#9FBBE0", label='VEE')
    ax1.fill_between(x_vals, y_ci, alpha=0.6, color="#FEDAB9", label='E_CI')
    ax1.plot(x_vals, y_vee, color="#9FBBE0", linewidth=2.5)
    ax1.plot(x_vals, y_ci, color="#FEDAB9", linewidth=2.5)
    ax1.set_xlabel("Energy (eV)", fontsize=50)
    ax1.set_ylabel("Probability Density", fontsize=50)
    ax1.tick_params(labelsize=30)
    ax1.set_xlim(0, 12)
    ax1.set_ylim(0, 1.05)
    ax1.legend(fontsize=40)
    plt.savefig(out_path / "CI_VEE_energy_kde.svg", transparent=True, bbox_inches='tight', facecolor='none', dpi=300)
    plt.close()

    # Plot delta E = VEE - CI_energy distribution
    fig = plt.figure(figsize=(18, 10*18/16))
    ax2 = fig.add_subplot(111)
    delta_e = VEE_energies - CI_energies
    x_min = -4
    x_max = 6
    kde_delta_e = stats.gaussian_kde(delta_e)
    x_delta_e = np.linspace(x_min, x_max, 400)
    y_delta_e = kde_delta_e(x_delta_e)
    y_delta_e = y_delta_e / y_delta_e.max()  # Normalize to [0, 1]
    ax2.fill_between(x_delta_e, y_delta_e, alpha=0.6, color="#df7373", label='VEE - E_CI')
    ax2.plot(x_delta_e, y_delta_e, color="#df7373", linewidth=2.5)
    ax2.set_xlabel("VEE - E_CI (eV)", fontsize=50)
    ax2.set_ylabel("Probability Density", fontsize=50)
    ax2.tick_params(labelsize=30)
    ax2.set_xlim(x_min, x_max)
    ax2.set_ylim(0, 1.05)
    ax2.legend(fontsize=40)
    plt.savefig(out_path / "delta_energy_kde.svg", transparent=True, bbox_inches='tight', facecolor='none', dpi=300)
    plt.close()

def plot_compound_type_energy(delta_energy_df, Compound_TYPES,out_path):
    """
    Draw energy distribution scatter plot for each compound type
    """

    for compound_type in Compound_TYPES:
        mol_list = Compound_TYPES[compound_type]
        if len(mol_list) < 2:
            continue
        label = compound_type.replace(' ','_')
        CI_energies = delta_energy_df.loc[delta_energy_df['Molecule'].isin(mol_list), 'CI_energy']
        VEE_energies = delta_energy_df.loc[delta_energy_df['Molecule'].isin(mol_list), 'VEE']

        # Plot VEE and CI energy
        fig = plt.figure(figsize=(18, 10*18/16))
        ax1 = fig.add_subplot(111)
        x_vals = np.linspace(0, 12, 500)
        kde_vee = stats.gaussian_kde(VEE_energies)
        kde_ci = stats.gaussian_kde(CI_energies)
        y_vee = kde_vee(x_vals)
        y_ci = kde_ci(x_vals)
        y_max = max(y_vee.max(), y_ci.max())
        y_vee = y_vee / y_max
        y_ci = y_ci / y_max
        ax1.fill_between(x_vals, y_vee, alpha=0.6, color="#9FBBE0", label='VEE')
        ax1.fill_between(x_vals, y_ci, alpha=0.6, color="#FEDAB9", label='E_CI')
        ax1.plot(x_vals, y_vee, color="#9FBBE0", linewidth=2.5)
        ax1.plot(x_vals, y_ci, color="#FEDAB9", linewidth=2.5)
        ax1.set_xlabel("Energy (eV)", fontsize=50)
        ax1.set_ylabel("Probability Density", fontsize=50)
        ax1.tick_params(labelsize=30)
        ax1.set_xlim(0, 12)
        ax1.set_ylim(0, 1.05)
        ax1.legend(fontsize=40)
        plt.savefig(out_path / f"CI_VEE_energy_kde_{label}.svg", transparent=True, bbox_inches='tight', facecolor='none', dpi=300)
        plt.close()

        # Plot delta E = VEE - CI_energy distribution
        fig = plt.figure(figsize=(18, 10*18/16))
        ax2 = fig.add_subplot(111)
        delta_e = VEE_energies - CI_energies
        x_min = -4
        x_max = 6
        kde_delta_e = stats.gaussian_kde(delta_e)
        x_delta_e = np.linspace(x_min, x_max, 500)
        y_delta_e = kde_delta_e(x_delta_e)
        y_delta_e = y_delta_e / y_delta_e.max()  # Normalize to [0, 1]
        ax2.fill_between(x_delta_e, y_delta_e, alpha=0.6, color="#df7373", label='VEE - E_CI')
        ax2.plot(x_delta_e, y_delta_e, color="#df7373", linewidth=2.5)
        ax2.set_xlabel("VEE - E_CI (eV)", fontsize=50)
        ax2.set_ylabel("Probability Density", fontsize=50)
        ax2.tick_params(labelsize=30)
        ax2.set_xlim(x_min, x_max)
        ax2.set_ylim(0, 1.05)
        ax2.legend(fontsize=40)
        plt.savefig(out_path / f"delta_energy_kde_{label}.svg", transparent=True, bbox_inches='tight', facecolor='none', dpi=300)
        plt.close()

def main():
    main_path = Path("real_work_dir")
    category_file = main_path/"final_all.csv"
    hdf5_path = main_path / "Final_property.hdf5"
    out_path = main_path / "ana_energy"
    out_path.mkdir(parents=True, exist_ok=True)
    energy_file = out_path / "CI_VEE_energies.csv"
    category_df = pd.read_csv(category_file,index_col=1)
    delta_energy_df = read_hdf5_file(hdf5_path, energy_file)
    Compound_TYPES = classify_compound_type(category_df,delta_energy_df['Molecule'].tolist())

    for compound_type in Compound_TYPES:
       print(f"{compound_type}: {len(Compound_TYPES[compound_type])} molecules")
    plot_all_energy(delta_energy_df, out_path)
    plot_compound_type_energy(delta_energy_df, Compound_TYPES,out_path)


if __name__ == "__main__":
    main()
