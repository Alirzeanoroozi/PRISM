import os
from utils import standard_data

def get_relative_asa_complex(template, save_directory):
    run_naccess(template, save_directory)
    
    relative_asa_complex = {}    
    with open(f"{save_directory}/{template}.rsa", 'r') as f:
        for rsaline in f.readlines():
            if rsaline.startswith("RES"):
                items = rsaline.split()
                residue_name = items[1]
                standard = standard_data(residue_name)
                if standard != -1:
                    chain = items[2]
                    residue_number = items[3]
                    absolute_asa = float(items[4])
                    relative_asa = absolute_asa * 100 / standard
                    relative_asa_complex[f"{residue_name}_{residue_number}_{chain}"] = relative_asa    
    
    return relative_asa_complex

def get_absolute_asa_complex(template, save_directory):
    run_naccess(template, save_directory)
    
    absolute_asa_complex = {}    
    with open(f"{save_directory}/{template}.asa", 'r') as f:
        for asaline in f.readlines():
            if asaline.startswith("ATOM"):
                items = asaline.split()
                residue_name = items[3]
                residue_number = items[4]
                chain = items[5]
                absolute_asa_complex[f"{residue_name}_{residue_number}_{chain}"] = [items[6], items[7], items[8]]
    
    return absolute_asa_complex

def run_naccess(template, save_directory):
    os.system(f"external_tools/naccess/naccess pdbs/{template[:4].lower()}.pdb")
    os.rename(f"{template[:4].lower()}.rsa", f"{save_directory}/{template}.rsa")
    os.rename(f"{template[:4].lower()}.asa", f"{save_directory}/{template}.asa")
    os.remove(f"{template[:4].lower()}.log")
