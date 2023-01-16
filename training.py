"""
    This script trains the objective function, using interatomic distance distributions that are computed from a dataset of known 3D structures (i.e. experimentally determined)

    To do so:
        1. It determines the number of models for each RNA PDB file
        2. Prepare the data for the calculation 
        3. Reports files preparation
        4. Distances computing
        5. Observed frequency calculation
        6. Reference frequency calculation
        7. Log ratio calculation
"""
import os
import pandas as pd
import math

from settings import pdb_lines, pdb_cols, base_pairs, base_list, col_names, intervals
from files_manager import *

def get_num_model(seq_ref, dir_path="PDB"):
    """Returns the number of models for any given RNA PDB file
    
    Parameters:
    seq_ref (str): The reference of the sequence to be checked
    dir_path (str): The directory path to be checked, default is "PDB"

    Returns:
    int: Returning the number of models
    """
    try:
        with open(f"{dir_path}/{seq_ref}.pdb", "r") as pdb_file:
            for line in pdb_file:
                if line[:6] == "NUMMDL":
                    num_model = int(line.split()[1])
                    return num_model

            return 1

    except:
        print(f"Can't open {dir_path}/{seq_ref}.pdb file")

def data_prep(seq_ref, dir_path="PDB", num_model=1):
    """Prepares the data source files for computing for any given RNA PDB file
    
    Parameters:
    seq_ref (str): The reference of the sequence to be checked
    dir_path (str): The directory path to be checked, default is "PDB"
    num_model (int)

    Returns:
    None: Generates .mdl files including needed data only
    """
    if num_model == 1:
        with open(f"pdb_models/{seq_ref}m{num_model}.mdl", "w") as mdl_file:
            with open(f"{dir_path}/{seq_ref}.pdb", "r") as pdb_file:
                for line in pdb_file:
                    if line.strip().split()[0] in pdb_lines and line.strip().split()[2] in pdb_cols:
                        line = line.replace("-", " -").strip().split()
                        line = ";".join(line) + "\n"
                        mdl_file.write(line)

    elif num_model > 1:
        # Get the starting line number of each model
        with open(f"{dir_path}/{seq_ref}.pdb", "r") as pdb_file:
            idx_list = []
            for idx, line in enumerate(pdb_file):
                if line.strip().split()[0] == "MODEL":
                    idx_list.append(idx)

                elif line.strip().split()[0] == "END":
                    idx_list.append(idx)
        
        # Generate mdl file for all models
        model = 1
        while model <= num_model:            
            with open(f"pdb_models/{seq_ref}m{model}.mdl", "w") as mdl_file:
                with open(f"{dir_path}/{seq_ref}.pdb", "r") as pdb_file:
                    for idx, line in enumerate(pdb_file):
                        if idx >= idx_list[model-1] and idx <= idx_list[model] and line.strip().split()[0] in pdb_lines and line.strip().split()[2] in pdb_cols:
                            line = line.replace("-", " -").strip().split()
                            line = ";".join(line) + "\n"
                            mdl_file.write(line)

            model += 1

def report_prep(rpt_dir="reports"):
    """Prepares the reports files by checking if a version is existing within the report directory or not.

    If reports exist, then nothing is done. However, if one of the four reports files is missing it will be created and initiated with 0 for all inputs
    
    Parameters:
    rpt_dir (str): The directory path where report will be saved, default is "reports"

    Returns:
    None: Check if reports exists, if not, they will be created with 0 as initial version
    """
    reports = ["distances", "obs_freq", "ref_freq", "log_ratio"]
    for report in reports:
        if os.path.exists(f"{rpt_dir}/{report}.txt"):
            print(f"The {report} report already exist.")

        else:
            col_names = ["Bases"]
            for interval in intervals:
                col_names.append(interval)

            with open(f"{rpt_dir}/{report}.txt", "w") as file_report:
                file_report.write(";".join(col_names) + "\n")

                for pair in base_list:
                    file_report.write(pair)

                    for interval in intervals:
                        file_report.write(f";{0}")
                    
                    file_report.write("\n")

    reports = ["tmp_dist"]
    for report in reports:
        if os.path.exists(f"{rpt_dir}/{report}.txt"):
            print(f"The {report} report already exist.")

        else:
            col_names = ["Bases"]
            for interval in intervals:
                col_names.append(interval)

            with open(f"{rpt_dir}/{report}.txt", "w") as file_report:
                file_report.write(";".join(col_names) + "\n")

                for pair in base_pairs:
                    file_report.write(pair)

                    for interval in intervals:
                        file_report.write(f";{0}")
                    
                    file_report.write("\n")

def calc_distances(seq_ref, dir_path="pdb_models", rpt_dir="reports"):
    """Calculate the distances and update the distance report file
    
    Parameters:
    seq_ref (str): The reference of the sequence to be checked
    dir_path (str): The directory path where the models are stored, default is "pdb_models"
    rpt_dir (str): The directory path where report will be saved, default is "reports"

    Returns:
    None: Update the temp. distances report file
    """
    distances_df = pd.read_csv(f"{rpt_dir}/tmp_dist.txt", sep=";")
    num_model = get_num_model(seq_ref)

    mdl_list = []

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

                        x_a = float(mdl_list[idx_1][6])
                        x_b = float(mdl_list[idx_2][6])
                        
                        y_a = float(mdl_list[idx_1][7])
                        y_b = float(mdl_list[idx_2][7])
                        
                        z_a = float(mdl_list[idx_1][8])
                        z_b = float(mdl_list[idx_2][8])

                        distance = math.sqrt(((x_a - x_b)**2) + ((y_a - y_b)**2)+ ((z_a - z_b)**2))

                        if distance >= 0. and distance <= 20.: 
                            dist_int = int(distance)
                        
                            if dist_int >= 0 and dist_int < 20:
                                col_idx = dist_int
                                col_name = intervals[col_idx]
                                # print(f"ROW {row_idx}, COL {col_idx}, PAIR {pair}, DIST {distance}")
                                distances_df.at[row_idx, col_name] += 1

                            elif dist_int == 20:
                                col_idx = 19
                                col_name = intervals[col_idx]
                                # print(f"ROW {row_idx}, COL {col_idx}, PAIR {pair}, DIST {distance}")
                                distances_df.at[row_idx, col_name] += 1

                        # print(f"{idx_1}, {mdl_list[idx_1]}")
                        # print(f"{idx_2}, {mdl_list[idx_2]}")
                        # print(f"{distance}")
    print(distances_df)
    distances_df.to_csv(f"{rpt_dir}/tmp_dist.txt", sep=";", index=False) 

def final_distance(rpt_dir="reports"):
    """Combine the similar base pairs (i.e. AU/UA, AG/GA..) within unique pair and generate the final distances report file
    
    Parameters:
    rpt_dir (str): The directory path where report will be saved, default is "reports"

    Returns:
    None: Update the temp. distances report file
    """
    try:
        with open(f"{rpt_dir}/distances.txt", "w") as final_dist:
            final_dist.write("Bases;0-1;1-2;2-3;3-4;4-5;5-6;6-7;7-8;8-9;9-10;10-11;11-12;12-13;13-14;14-15;15-16;16-17;17-18;18-19;19-20\n")

            for pair in base_list:
                line_list = []

                with open(f"{rpt_dir}/tmp_dist.txt") as tmp_dist:
                    for tmp_line in tmp_dist:
                        if tmp_line.strip().split(";")[0] == pair or tmp_line.strip().split(";")[0] == pair[::-1]:
                            line_list.append(tmp_line.strip().split(";"))

                    if len(line_list) == 1:
                        # print(";".join(line_list[0]))
                        final_dist.write(";".join(line_list[0]) + "\n")

                    if len(line_list) == 2:
                        # print(pair, end=";")
                        final_dist.write(pair + ";")

                        for interval in intervals:
                            idx = intervals.index(interval) + 1
                            x = int(line_list[0][idx]) + int(line_list[1][idx])

                            if interval == "19-20":
                                # print(f"{x}", end="")
                                final_dist.write(f"{x}")
                            else:
                                # print(f"{x}", end=";")
                                final_dist.write(f"{x}" + ";")

                        # print()
                        final_dist.write("\n")

                        # print(line_list[0], '\n', line_list[1])

    except:
        print("Something wrong with the distances report, please try to re-run the code from the begging!") 

def calc_obs_freq(rpt_dir="reports"):
    """Calculate the observed frequency and update the report file
    
    Parameters:
    rpt_dir (str): The directory path where report will be saved, default is "reports"

    Returns:
    None: Update the observed frequency report file
    """
    try:
        distances_df = pd.read_csv(f"{rpt_dir}/distances.txt", sep=";")
        distances_df.loc[:,"Row_Total"] = distances_df.sum(numeric_only=True, axis=1)

        obs_freq_df = pd.read_csv(f"{rpt_dir}/obs_freq.txt", sep=";")

        for row_idx, pair in enumerate(base_list):
            for interval in intervals:
                obs_freq_df.at[row_idx, interval] = distances_df.at[row_idx, interval] / distances_df.at[row_idx, "Row_Total"]

        print(obs_freq_df)
        obs_freq_df.to_csv(f"{rpt_dir}/obs_freq.txt", sep=";", index=False)

    except:
        print("Something wrong with the distances report, please try to re-run the code from the begging!")

def calc_ref_freq(rpt_dir="reports"):
    """Calculate the reference frequency and update the report file
    
    Parameters:
    rpt_dir (str): The directory path where report will be saved, default is "reports"

    Returns:
    None: Update the reference frequency report file
    """
    try:
        distances_df = pd.read_csv(f"{rpt_dir}/distances.txt", sep=";")
        distances_df.loc["Column_Total"]= distances_df.sum(numeric_only=True, axis=0)

        ref_freq_df = pd.read_csv(f"{rpt_dir}/ref_freq.txt", sep=";")

        for row_idx, pair in enumerate(base_list):
            for interval in intervals:
                ref_freq_df.at[row_idx, interval] = distances_df.at[row_idx, interval] / distances_df.at["Column_Total", interval]

        print(ref_freq_df)
        ref_freq_df.to_csv(f"{rpt_dir}/ref_freq.txt", sep=";", index=False)

    except:
        print("Something wrong with the distances report, please try to re-run the code from the begging!")

def calc_log_ratio(rpt_dir="reports"):
    """Calculate the log ratio and update the report file
    
    Parameters:
    rpt_dir (str): The directory path where report will be saved, default is "reports"

    Returns:
    None: Update the log ratio report file
    """
    try:
        log_ratio_df = pd.read_csv(f"{rpt_dir}/log_ratio.txt", sep=";")

        ref_freq_df = pd.read_csv(f"{rpt_dir}/ref_freq.txt", sep=";")
        obs_freq_df = pd.read_csv(f"{rpt_dir}/obs_freq.txt", sep=";")

        for row_idx, pair in enumerate(base_list):
            for interval in intervals:
                log_ratio_df.at[row_idx, interval] = -1 * math.log10(obs_freq_df.at[row_idx, interval] / ref_freq_df.at[row_idx, interval])

        log_ratio_df = log_ratio_df.fillna(10)
        print(log_ratio_df)
        log_ratio_df.to_csv(f"{rpt_dir}/log_ratio.txt", sep=";", index=False)

    except:
        print("Something wrong with the frequencies reports, please try to re-run the code from the begging!")

def training_run(seq_ref, dir_path="pdb_models"):
    """The main training script, it trains the objective function, using interatomic distance distributions that are computed from a dataset of known 3D structures (i.e. experimentally determined)
    
    Parameters:
    seq_ref (str): The reference of the sequence to be checked
    dir_path (str): The directory path to store the models, default is "pdb_models"

    Returns:
    None: Launch the training script to perform the requested computing
    """
    print(f"Training using {seq_ref} started")

    num_model = get_num_model(seq_ref)
    data_prep(seq_ref, num_model=num_model)

    print("Report files preparation...")
    report_prep()

    print("Distances calculation...")
    calc_distances(seq_ref)
    final_distance()

    print("Observed frequency calculation...")
    calc_obs_freq()

    print("Reference frequency calculation...")
    calc_ref_freq()

    print("Log ratio calculation...")
    calc_log_ratio()

if __name__ == "__main__":
    print("Welcome to the Training Script...")