"""
    This script generates the different plots. it will plot the interaction profiles using the matplotlib library of the score as a function of the distance.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from settings import pdb_lines, pdb_cols, base_pairs, base_list, col_names, intervals
from files_manager import *

def plot(rpt_dir="reports", plt_dir="plot"):
    """Plot the interaction profiles: the score as a function of the distance.
    Matplotlib library was used to generate the different graphs, which are saved to 
    
    Parameters:
    rpt_dir (str): The directory path where reports are saved, default is "reports"
    plt_dir (str): The directory where the plot images will be saved

    Returns:
    None: Generate and save the interaction profiles
    """

    log_ratio_df = pd.read_csv(f"{rpt_dir}/log_ratio.txt", sep=";")
    dfT = log_ratio_df.set_index("Bases").T

    for pair in base_list:
        plt.rc('grid', linestyle=':')
        dfT.plot(y = pair)

        # Add title and axis names
        plt.title(f"Interaction Profiles Plot Results of {pair} Base-Pair")
        plt.xlabel(r"Distances [$\AA$]")
        plt.ylabel("Pseudo-Energy")

        # Create names on the x axis
        plt.xticks(range(len(intervals) + 1), range(len(intervals) + 1), rotation=45)

        # Add X-Axis Grid
        plt.grid(axis='x')

        # Save the plot to a PNG file
        plt.savefig(f"{plt_dir}/{pair}.png")
    
    # Show the plot
    plt.show()



if __name__ == "__main__":
    print("Welcome to the Plot Script...")