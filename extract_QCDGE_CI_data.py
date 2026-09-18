import h5py
from pathlib import Path
import numpy as np
import pandas as pd
from tqdm import tqdm
from ase.data import chemical_symbols
import argparse

def extract_data_from_hdf5_file(hdf5_file,subgroup,energy_file,coord_file):
    """
    Extract data(energy or coordinates) from hdf5 file and save to csv or xyz file.
    """
    property_type = ["energy", "coordinates"]
    property_choose = int(input("Choose 0 for energy, 1 for coordinates: "))
    with h5py.File(hdf5_file, 'r') as hdf:
        groups =list(hdf.keys())
        for group in tqdm(groups, desc="Processing groups"):
            group_data = hdf[group]
            labels = group_data.attrs["labels"]
            atoms = [chemical_symbols[atomic_number] for atomic_number in labels]
            subgroup_data = group_data[subgroup]
            property = property_type[property_choose]
            data = subgroup_data[property][()]
            if isinstance(data, bytes):
                data = data.decode('utf-8')
            if property == "energy":
                energy = data
                save_energy_to_csv(group,subgroup,energy,energy_file)
            elif property == "coordinates":
                coordinates = data
                save_data_to_xyz(coordinates, atoms, group, coord_file)

def save_data_to_xyz(coordinates, atoms,group_name,filename):
        num_atoms = len(coordinates)
        with open(filename, 'a') as f:
            f.write(f"{num_atoms}\n")
            f.write(f"Group: {group_name}\n")
            for i in range(num_atoms):
                f.write(f"{atoms[i]} {coordinates[i, 0]} {coordinates[i, 1]} {coordinates[i, 2]}\n")

def save_energy_to_csv(group,subgroup,energy,filename):
    with open(filename, 'a') as f:
        f.write(f"{group},{subgroup},{energy}\n")


def extract_data_from_QCDGE_CI():
    parser = argparse.ArgumentParser(description='Extract data from hdf5 file')
    parser.add_argument('--hdf5_file', type=Path, help='Path to the hdf5 file')
    parser.add_argument('--main_path', type=Path, help='Path to the main directory')
    args = parser.parse_args()
    hdf5_file = args.hdf5_file
    main_path = args.main_path
    subgroups=["Ground", "CI"]
    choose = int(input("Choose 0 for Ground state properties, 1 for CI state properties: "))
    energy_file = main_path / f"energy_{subgroups[choose]}.csv"
    coord_file = main_path / f"coords_{subgroups[choose]}.xyz"
    extract_data_from_hdf5_file(hdf5_file,subgroups[choose],energy_file,coord_file)

if __name__ == "__main__":
    extract_data_from_QCDGE_CI()

    #command line
    #python extract_QCDGE_CI_data.py --hdf5_file input_hdf5_file_path --main_path save_dir_path