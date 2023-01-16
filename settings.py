# Code details
__author__ = "Assia B."
__version__ = "1.0.0"

# Default directories
dir_list = ["PDB", "pdb_models", "plot", "reports"]

# Lines to keep pn pdb files
pdb_lines = ["ATOM"]
pdb_cols = ["C3'"]

# Reports initial settings
base_pairs = ["AA" ,"AU", "AC", "AG", "UU", "UC", "UG", "CC","CG", "GG", "GC","GU","CU","GA","CA","UA"]
base_list = ["AA", "AU", "AC", "AG", "UU", "UC", "UG", "CC", "CG", "GG"]

intervals = ["0-1", "1-2", "2-3", "3-4", "4-5", "5-6", "6-7", "7-8", "8-9", "9-10", "10-11", "11-12", "12-13", "13-14", "14-15", "15-16", "16-17", "17-18", "18-19", "19-20"]

# DataFrame column names
col_names = ["record_type", "atom_num", "atom", "base", "chain_id", "residue_num", "coor_x", "coor_y", "coor_z", "occupancy", "temp_factor", "element_name"]