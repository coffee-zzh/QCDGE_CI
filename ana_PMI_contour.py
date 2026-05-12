from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem import Descriptors3D
from pathlib import Path
from typing import List, Tuple, Dict, Any, Optional
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
from scipy.stats import gaussian_kde
from matplotlib.path import Path as mlPath
import h5py
from ase.data import chemical_symbols
import time


# Constants
ATOMIC_MASSES = {'H': 1.007825, 'C': 12.011, 'N': 14.0067, 'O': 15.9994, 'F': 18.9984032}
SHAPE_THRESHOLD = 0.05
IDEAL_SHAPES = {
    'rod': (0.0, 1.0),
    'disk': (0.5, 0.5),
    'sphere': (1.0, 1.0)
}
TRIANGLE_VERTICES = np.array([[0, 1], [1, 1], [0.5, 0.5], [0, 1]])


def get_masses(symbols: List[str]) -> List[float]:
    """
    Get masses of atoms in a molecule.

    Args:
        symbols: List of atomic symbols

    Returns:
        List of atomic masses corresponding to the symbols
    """
    masses = []
    for symbol in symbols:
        if symbol not in ATOMIC_MASSES:
            raise ValueError(f"Unknown atomic symbol: {symbol}")
        masses.append(ATOMIC_MASSES[symbol])
    return masses


def cal_xyz_Pmi(coord: np.ndarray, masses: List[float]) -> Tuple[float, float, float, float, float]:
    """
    Calculate PMI (Principal Moments of Inertia) and NPR (Normalized Principal Ratios)
    for a molecule by coordinates and masses.

    Args:
        coord: Array of shape (N, 3) containing atomic coordinates
        masses: List of atomic masses

    Returns:
        Tuple of (pmi1, pmi2, pmi3, npr1, npr2)
        where pmi1 <= pmi2 <= pmi3 are the sorted principal moments
        and npr1 = pmi1/pmi3, npr2 = pmi2/pmi3

    Raises:
        ValueError: If input dimensions are incorrect or total mass is zero
    """
    coord = np.asarray(coord)
    masses = np.asarray(masses)

    if coord.ndim != 2 or coord.shape[1] != 3:
        raise ValueError(f"coord must be shape (N, 3), got {coord.shape}")
    if len(masses) != len(coord):
        raise ValueError(f"masses length ({len(masses)}) must match coord length ({len(coord)})")

    total_mass = np.sum(masses)
    if total_mass <= 0:
        raise ValueError("Total mass must be positive")

    # Center of mass
    center_of_mass = np.sum(coord * masses[:, np.newaxis], axis=0) / total_mass
    coord_com = coord - center_of_mass

    # Inertia tensor components
    Ixx = np.sum(masses * (coord_com[:, 1]**2 + coord_com[:, 2]**2))
    Iyy = np.sum(masses * (coord_com[:, 0]**2 + coord_com[:, 2]**2))
    Izz = np.sum(masses * (coord_com[:, 0]**2 + coord_com[:, 1]**2))
    Ixy = -np.sum(masses * coord_com[:, 0] * coord_com[:, 1])
    Ixz = -np.sum(masses * coord_com[:, 0] * coord_com[:, 2])
    Iyz = -np.sum(masses * coord_com[:, 1] * coord_com[:, 2])

    inertia_tensor = np.array([[Ixx, Ixy, Ixz],
                               [Ixy, Iyy, Iyz],
                               [Ixz, Iyz, Izz]])

    eigvals = np.linalg.eigvalsh(inertia_tensor)
    pmi1, pmi2, pmi3 = np.sort(eigvals)

    # Avoid division by zero
    if pmi3 <= 0:
        raise ValueError("Largest principal moment must be positive")

    npr1 = pmi1 / pmi3
    npr2 = pmi2 / pmi3

    return pmi1, pmi2, pmi3, npr1, npr2


def classify_shape(npr1: float, npr2: float) -> str:
    """
    Classify shape of a molecule based on NPR values.

    Args:
        npr1: Normalized principal ratio I1/I3
        npr2: Normalized principal ratio I2/I3

    Returns:
        Shape classification: 'rod', 'disk', 'sphere', or 'other'
    """
    distances = {}
    for shape, (ideal_npr1, ideal_npr2) in IDEAL_SHAPES.items():
        dist = np.sqrt((npr1 - ideal_npr1)**2 + (npr2 - ideal_npr2)**2)
        distances[shape] = dist

    min_shape = min(distances, key=distances.get)
    min_distance = distances[min_shape]

    if min_distance < SHAPE_THRESHOLD:
        return min_shape
    else:
        return 'other'


def _draw_triangle(ax: plt.Axes) -> None:
    """
    Draw the PMI triangle boundary on a matplotlib Axes.

    Args:
        ax: Matplotlib Axes object to draw on
    """
    ax.fill(TRIANGLE_VERTICES[:, 0], TRIANGLE_VERTICES[:, 1], alpha=0.1, color='lightblue')
    ax.plot(TRIANGLE_VERTICES[:, 0], TRIANGLE_VERTICES[:, 1], 'k-', linewidth=2)


def cal_PMI(hdf5_path: Path, label: str) -> Tuple[List[List[float]], List[List[float]], List[str], List[str], List[str], List[str], List[str]]:
    """
    Calculate PMI and NPR for a list of S0-S1 CI molecules from HDF5 file.

    Args:
        hdf5_path: Path to HDF5 file containing CI data

    Returns:
        Tuple of (pmi_list, npr_list, disk_list, rod_list, sphere_list, other_list)
    """
    pmi_list = []
    npr_list = []
    disk_list = []
    rod_list = []
    sphere_list = []
    other_list = []
    error_list = []

    with h5py.File(hdf5_path, 'r') as f:
        groups = list(f.keys())
        

        for group in tqdm(groups, desc=f"Processing {label} molecules"):
            try:
                group_data = f[group]
                labels = group_data.attrs['labels'][()]
                atoms = [chemical_symbols[atomic_number] for atomic_number in labels]
                coord_data = group_data[label]
                coord = coord_data['coordinates'][()]
                masses = get_masses(atoms)

                pmi1, pmi2, pmi3, npr1, npr2 = cal_xyz_Pmi(coord, masses)
                pmi_list.append([pmi1, pmi2, pmi3])
                npr_list.append([npr1, npr2])

                shape = classify_shape(npr1, npr2)
                if shape == 'rod':
                    rod_list.append(group)
                elif shape == 'disk':
                    disk_list.append(group)
                elif shape == 'sphere':
                    sphere_list.append(group)
                else:
                    other_list.append(group)

            except Exception as e:
                print(f"Error processing group {group}: {e}")
                error_list.append(group)
                continue

    return pmi_list, npr_list, disk_list, rod_list, sphere_list, other_list, error_list


def save_PMI_results(
    out_path: Path,
    npr_list: List[List[float]],
    pmi_list: List[List[float]],
    disk_list: List[str],
    rod_list: List[str],
    sphere_list: List[str],
    other_list: List[str],
    label: str
    ) -> None:
    """
    Save CI PMI analysis results to CSV files.

    Args:
        out_path: Output directory path
        npr_list: List of [npr1, npr2] pairs
        pmi_list: List of [pmi1, pmi2, pmi3] triplets
        disk_list: List of group IDs classified as disk
        rod_list: List of group IDs classified as rod
        sphere_list: List of group IDs classified as sphere
        other_list: List of group IDs classified as other
    """
    # Build shape dictionary for O(1) lookups
    shape_map = {}
    for group in disk_list:
        shape_map[group] = 'disk'
    for group in rod_list:
        shape_map[group] = 'rod'
    for group in sphere_list:
        shape_map[group] = 'sphere'
    for group in other_list:
        shape_map[group] = 'other'

    all_groups = disk_list + rod_list + sphere_list + other_list
    data = []

    for i, group in enumerate(all_groups):
        if i < len(npr_list) and i < len(pmi_list):
            nprs = npr_list[i]
            pmi = pmi_list[i]
            shape = shape_map.get(group, 'unknown')

            data.append({
                'Group': group,
                'NPR1': nprs[0],
                'NPR2': nprs[1],
                'PMI1': pmi[0],
                'PMI2': pmi[1],
                'PMI3': pmi[2],
                'Shape': shape
            })

    df_results = pd.DataFrame(data)

    csv_path = out_path / f'{label}_shape_analysis_results.csv'
    df_results.to_csv(csv_path, index=False)
    print(f"{label} shape analysis results saved to: {csv_path}")



def plot_pmi_figure(npr_list: List[List[float]], out_path: Path) -> None:
    """
    Plot PMI contour figure with KDE density.

    Args:
        npr_list: List of [npr1, npr2] pairs
        out_path: Output directory path
    """
    plt.figure(figsize=(18, 10 * 18 / 16))
    ax = plt.gca()

    # Draw triangle boundary
    _draw_triangle(ax)

    if npr_list:
        data = np.array(npr_list)
        x = data[:, 0]
        y = data[:, 1]

        xi = np.linspace(-0.05, 1.05, 600)
        yi = np.linspace(0.45, 1.05, 600)
        Xi, Yi = np.meshgrid(xi, yi)

        xy = np.vstack([x, y])
        kde = gaussian_kde(xy)

        positions = np.vstack([Xi.ravel(), Yi.ravel()])
        Zi = kde(positions).reshape(Xi.shape)

        triangle_path = mlPath(TRIANGLE_VERTICES)

        points = np.vstack([Xi.ravel(), Yi.ravel()]).T
        mask = ~triangle_path.contains_points(points).reshape(Xi.shape)
        Zi_masked = Zi.copy()
        Zi_masked[mask] = np.nan

        contour = plt.contour(Xi, Yi, Zi_masked,
                                levels=np.linspace(0, 15, 20),
                                cmap='viridis_r',
                                alpha=0.7,
                                vmin= 0,
                                vmax= 15,
                                extend='both',
                                linewidths=2)
   

        cbar = plt.colorbar(contour)
        cbar.ax.tick_params(labelsize=40)
        cbar.set_label('Probability Density', fontsize=50)
        

    plt.xlabel('I₁/I₃', fontsize=50)
    plt.ylabel('I₂/I₃', fontsize=50)
    plt.title(f'PMI Analysis ', fontsize=50)
    plt.xlim(-0.05, 1.05)
    plt.ylim(0.45, 1.05)
    plt.tick_params(labelsize=30)

    plt.tight_layout()
    plt.savefig(out_path / 'pmi_density_contour_line.svg', dpi=300, bbox_inches='tight')
    plt.close()



def analyze_Ground_PMI(hdf5_path: Path, out_path: Path) -> None:
    """
    Analyze ground state PMI from SMILES/InChI CSV file.

    Args:
        smiles_path: Path to CSV file with Smiles_rdkit and Inchi_rdkit columns
        out_path: Output directory path
    """
    pmi_list, npr_list, disk_list, rod_list, sphere_list, other_list, error_list= cal_PMI(hdf5_path,'Ground')
    plot_pmi_figure(npr_list, out_path)
    save_PMI_results(out_path, npr_list, pmi_list, disk_list, rod_list, sphere_list, other_list,label='Ground')
    if error_list:
        error_df = pd.DataFrame(error_list, columns=['Group'])
        error_df.to_csv(out_path / 'error_groups.csv', index=False)
        print(f"Total {len(error_list)} groups failed to process.")


def analyze_CI_PMI(hdf5_path: Path, out_path: Path) -> None:
    """
    Analyze CI state PMI from HDF5 file.

    Args:
        hdf5_path: Path to HDF5 file with CI data
        out_path: Output directory path
    """
    pmi_list, npr_list, disk_list, rod_list, sphere_list, other_list,error_list = cal_PMI(hdf5_path,'CI')
    plot_pmi_figure(npr_list, out_path)
    save_PMI_results(out_path, npr_list, pmi_list, disk_list, rod_list, sphere_list, other_list,label='CI')
    if error_list:
        error_df = pd.DataFrame(error_list, columns=['Group'])
        error_df.to_csv(out_path / 'error_groups.csv', index=False)
        print(f"Total {len(error_list)} groups failed to process.")



def main() -> None:
    """Main function to run PMI analysis pipeline."""
    main_path = Path('real_work_dir')
    hdf5_path = main_path / 'Final_property.hdf5'
    out_path = main_path / 'ana_PMI'
    out_path.mkdir(parents=True, exist_ok=True)
    start_time = time.time()
    # Ground state analysis
    ground_pmi_path = out_path / 'Ground_PMI'
    ground_pmi_path.mkdir(parents=True, exist_ok=True)
    analyze_Ground_PMI(hdf5_path, ground_pmi_path)

    # CI state analysis
    ci_pmi_path = out_path / 'CI_PMI'
    ci_pmi_path.mkdir(parents=True, exist_ok=True)
    analyze_CI_PMI(hdf5_path, ci_pmi_path)
    end_time = time.time()
    print(f"Total time: {end_time - start_time:.2f} seconds")



if __name__ == '__main__':
    main()
