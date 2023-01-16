"""
    This code was developed for the RNA folding energy estimation
    Three main steps were identified: (i) training, (ii) plot and (iii) scoring

    This version was developed for command line user
"""
from sys import exit
import argparse

from settings import __version__, __author__

from files_manager import *
from training import *
from plot import *
from evaluation import *

parser = argparse.ArgumentParser(description=f"Welcome to the RNA Folding Energy Estimator {__version__}. Created by {__author__}")

parser.add_argument('--seq', type=str, help="Run the code on specific RNA using its sequence reference")
parser.add_argument('-s', type=str, help="Run the code on specific RNA using its sequence reference")

parser.add_argument('--list', type=str, help="Run the code on list of RNA where their sequence references are stored on a specific file")
parser.add_argument('-l', type=str, help="Run the code on list of RNA where their sequence references are stored on a specific file")

args = parser.parse_args()

def main():
    dir_prep()

    if args.seq:
        seq_ref = args.seq.upper()
        get_pdb(seq_ref)
    
    elif args.s:
        seq_ref = args.s.upper()
        get_pdb(seq_ref)

    elif args.list:
        seq_list = args.list
        
        try:
            with open(seq_list) as my_list:
                pdb_list = my_list.readlines()

                for pdb in pdb_list:
                    get_pdb(pdb.strip())
        
        except:
            print(f"The file {seq_list} does not exist")

    elif args.l:
        seq_list = args.l
        
        try:
            with open(seq_list) as my_list:
                pdb_list = my_list.readlines()

                for pdb in pdb_list:
                    get_pdb(pdb.strip())
        
        except:
            print(f"The file {seq_list} does not exist")
    
    else:
        print("Please select one of the available options. For more details please check the help.")

        exit
    
    # Generating a list of all structures and another one for RNA only
    pdb_list = get_pdb_list()
    rna_list = []
    for pdb in pdb_list:
        if is_rna(pdb):
            rna_list.append(pdb)
            
    print(f"The following structures are available for analysis: {pdb_list}")
    print(f"The RNA ones (which will be used) are: {rna_list}")

    # Run the training script for each RNA
    try:
        for rna in rna_list:
            training_run(rna)

        plot()
    
    except:
        print("No RNA file for the analysis")

if __name__ == "__main__":
    main()