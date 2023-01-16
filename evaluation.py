"""
    This script is partially similar to the training one, as it will compute all the distances for a given structure (same thresholds: 20 A and i, i+4). For each distance, a scoring value will be computed, using a linear interpolation. By summing all these scores, the script will calculate the estimated Gibbs free energy of the evaluated RNA conformation.
"""
import os
import pandas as pd
import math

from settings import pdb_lines, pdb_cols, base_pairs, base_list, col_names, intervals

from files_manager import *
from training import *

def linear_interpolation(seq_ref, dir_path="pdb_models", rpt_dir="reports"):
    """Calculate of the Gibbs energy based on the distances and return the calculated value
    
    Parameters:
    seq_ref (str): The reference of the sequence to be checked
    dir_path (str): The directory path where the models are stored, default is "pdb_models"
    rpt_dir (str): The directory path where report will be saved, default is "reports"

    Returns:
    gibbs_energy (float): Returns the gibbs energy
    """
    train_score_df = pd.read_csv(f"{rpt_dir}/log_ratio.txt", sep=";")
    num_model = get_num_model(seq_ref)

    mdl_list = []
    gibbs_energy = 0

    for idx in range(num_model):
        print(f"Working on Seq. {seq_ref} - Model No. {idx + 1}")
        with open(f"{dir_path}/{seq_ref}m{idx + 1}.mdl", "r") as mdl_file:
            for line in mdl_file:
                mdl_list.append(line.strip().split(";"))
            
        for idx_1 in range(0, len(mdl_list) - 4):
            for idx_2 in range(idx_1 + 4, len(mdl_list)):
                # print(f"{idx_1} - {idx_2}")
                if mdl_list[idx_2][4] == mdl_list[idx_1][4]:
                    pair = f"{mdl_list[idx_1][3]}{mdl_list[idx_2][3]}"
                    # print(pair)

                    if pair in base_pairs:
                        row_idx = base_pairs.index(pair)
                        try:
                            row_idx = base_list.index(pair)
                        except:
                            row_idx = base_list.index(pair[::-1])
                        
                        x_a = float(mdl_list[idx_1][6])
                        x_b = float(mdl_list[idx_2][6])
                        
                        y_a = float(mdl_list[idx_1][7])
                        y_b = float(mdl_list[idx_2][7])
                        
                        z_a = float(mdl_list[idx_1][8])
                        z_b = float(mdl_list[idx_2][8])

                        distance = math.sqrt(((x_a - x_b)**2) + ((y_a - y_b)**2)+ ((z_a - z_b)**2))

                        x_1 = math.floor(distance)
                        x_2 = math.ceil(distance)
                        
                        if distance >= 0. and distance <= 20.:
                            dist_int = int(distance)
                            col_idx = dist_int
                            col_name = intervals[col_idx]                    
                            
                            y_1 = train_score_df.iloc[row_idx, col_idx]
                            y_2 = train_score_df.iloc[row_idx, col_idx]
                           
                            energy = y_1 + (distance - x_1) * (y_2 - y_1) / (x_2 - x_1)
                            dist_int = int(distance)

                            gibbs_energy += energy
                        
    print(f"The Gibbs energy is: {gibbs_energy}")
    return gibbs_energy

def evaluation_run(seq_ref):
    """The main evaluation script, it computes the scoring value using a linear interpolation
    
    Parameters:
    seq_ref (str): The reference of the sequence to be checked

    Returns:
    None: Launch the training script to perform the requested computing
    """
    print(f"Evaluation script for {seq_ref} started")

    get_pdb(seq_ref)

    if is_rna(seq_ref):
        
        num_model = get_num_model(seq_ref)
        data_prep(seq_ref, num_model=num_model)

        print("Report files preparation...")
        report_prep()

        print("Distances calculation...")
        calc_distances(seq_ref)
        final_distance()

        print("Score calculation...")
        linear_interpolation(seq_ref)
        
    
    else:
        print(f"The {seq_ref} PDB file is not an RNA.")
    

if __name__ == "__main__":
    print("Welcome to the Evaluation Script...")

    #dir_prep()
    evaluation_run("4P5J")