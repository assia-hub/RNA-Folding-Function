"""
    This Files Manager contains all needed functions allowing the assets, directories and files management tools
    
    1. Clean all files within specific directory
    2. Download PDB automatically using the sequence reference
    3. Check if the PDB file is for an RNA 
    4. Get the list of available PDB files for the analysis with the PDB directory
    5. Copy PDB files from a source directory to a destination one
"""
import os
import shutil
import wget

from settings import dir_list

def dir_clean(dir_path):
    """Remove all existing files inside the given directory
    
    Parameters:
    path (str): The path of the directory to be cleaned

    Returns:
    None
    """
    for file in os.listdir(dir_path):
        os.remove(os.path.join(dir_path, file))
        print(f"The {file} file was successfully removed.")

def get_pdb(seq_ref, dir_path="PDB"):
    """Downloads PDB file using the sequence reference if exists
    
    Parameters:
    seq_ref (str): The reference of the sequence to be downloaded
    dir_path (str): The directory path where files will be downloaded, default is "PDB"

    Returns:
    None: Downloads the PDB file to PDB directory
    """
    seq_ref = seq_ref.upper()

    try:
        url = f"https://files.rcsb.org/view/{seq_ref}.pdb"
        response = wget.download(url, f"{dir_path}/{seq_ref}.pdb")

        print("... Download done!")

    except:
        print(f"The given {seq_ref} sequence ID does not exit!")

def is_rna(seq_ref, dir_path="PDB"):
    """Check if the structure is RNA or not, and returns TRUE or FALSE
    
    Parameters:
    seq_ref (str): The reference of the sequence to be checked
    dir_path (str): The directory path to be checked, default is "PDB"

    Returns:
    bool: Returning True if RNA and False if not
    """
    try:
        with open(f"{dir_path}/{seq_ref}.pdb") as pdb_file:
            # file_header = pdb_file.readline().strip().replace("/", " ").replace("-", " ").split()
            file_header = pdb_file.readline().split()

            if "RNA" in file_header:
                return True
            else:
                return False
    
    except:
        print(f"Can not open the {dir_path}/{seq_ref}.pdb file")
        return False

def get_pdb_list(dir_path="PDB"):
    """Returns the complete list of the pdb files within the PDB directory available for the analysis
    
    Parameters:
    dir_path (str): The directory path to be checked, default is "PDB"

    Returns:
    list: Returning list of strings with pdb file names
    """
    pdb_list = []

    try:
        for file in os.listdir(dir_path):
            if os.path.isfile(os.path.join(dir_path, file)):
                if file.endswith(".pdb") or file.endswith(".PDB"):
                    pdb_list.append(file.replace(".pdb", "").replace(".PDB", ""))

    except:
        print(f"The directory {dir_path} does not exist.")

    return pdb_list

def dir_prep():
    """Create and clean needed directories

    PDB: To include pdb files which are used for the training script
    
    Parameters:
    None

    Returns:
    None
    """

    for dir in dir_list:
        if os.path.exists(dir):
            shutil.rmtree(dir)
        
        os.mkdir(dir)

def cp_pdb(src_dir, dir_path="PDB"):
    """Copy PDB files from the src_dir to dir_path
    
    Parameters:
    src_dir (str): The directory path including the PDB files
    dir_path (str): The directory path where PDB files will be copied, default is "PDB"

    Returns:
    None: Copy the PDB files to PDB directory
    """
    try:
        for file in os.listdir(src_dir):
            if file.endswith(".pdb") or file.endswith(".PDB"):
                shutil.copy(f"{src_dir}/{file}", f"{dir_path}/{file}")

    except:
        print(f"The directory {src_dir} does not exist.")

if __name__ == "__main__":
    print("Welcome to the File Manager...")