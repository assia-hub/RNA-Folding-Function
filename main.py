"""
    This code was developed for the RNA folding energy estimation
    Three main steps were identified: (i) training, (ii) plot and (iii) scoring
"""
from sys import exit

from settings import __version__, __author__

from files_manager import *
from training import *
from plot import *
from evaluation import *

def main():
    """The main function which manage the RNA Folding project scripts

    Parameters: None

    Returns: None
    """

    print(f"Welcome to the RNA Folding Energy Estimator {__version__}. Created by {__author__}\n")

    # STEP 0. Directories Preparation
    print("Before starting, do you want to clean previous data?")
    while True:
        user_choice = input("Please enter [1] for YES or [2] for NO: ")

        if user_choice == "1":
            dir_prep()
            break

        elif user_choice == "2":
            break

        else:
            print(f"\nYou entered [{user_choice}], this choice is not allowed.")
            print("Please select from the following list.\n")  

    # STEP 1. Call the Training Script Based on the User Selection
    print("Please indicate how you would like to select the PDB files?")
    
    while True:
        print("The PDF files can be selected by: ")
        print("[1] - Indicating the sequence reference")
        print("[2] - Indicating a file containing a list of sequences references")
        print("[3] - Indicating a directory containing PDB files")
        print("[4] - Just continue using previous files")
        print("[5] - Stop the program and quit")

        user_choice = input("\nPlease enter your choice and press enter: ")

        if user_choice == "1":
            seq_ref = input("Please enter the SEQ REF and press enter (e.g. 1ATW): ")
            get_pdb(seq_ref)

            break

        elif user_choice == "2":
            seq_list = input("Please enter the path to the file including the list of PDB files to be downloaded or just press enter for the default file (pdb_list.txt): ")

            if seq_list == "":
                seq_list = "pdb_list.txt"

            with open(seq_list) as my_list:
                pdb_list = my_list.readlines()

                for pdb in pdb_list:
                    get_pdb(pdb.strip())

            break

        elif user_choice == "3":
            seq_dir = input("Please enter the path to the directory including the PDB files to be copied or just press enter to use the default directory (pdb_files): ")

            if seq_dir == "":
                seq_dir = "pdb_files"

            cp_pdb(seq_dir)

            break

        elif user_choice == "4":
            break 

        elif user_choice == "5":
            exit()

        else:
            print(f"\nYou entered [{user_choice}], this choice is not allowed.")
            print("Please select from the following list.\n")  
    
    # Generating a list of all structures and another one for RNA only
    pdb_list = get_pdb_list()
    rna_list = []
    for pdb in pdb_list:
        if is_rna(pdb):
            rna_list.append(pdb)
            
    print(f"The following structures are available for analysis: {pdb_list}")
    print(f"The RNA ones (which will be used) are: {rna_list}")

    # Run the training script for each RNA
    for rna in rna_list:
        training_run(rna)

    plot()

if __name__ == "__main__":
    main()