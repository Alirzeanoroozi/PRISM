import os
import pickle
from Bio.PDB import PDBParser

os.makedirs("alignment", exist_ok=True)

def align(queries, templates):
    for protein in queries:
        for template in templates:
            for chain in template[4:]:
                protein_path = f"surface_extraction/{protein}.asa.pdb"
                interface_path = f"template/interface/{template}_{chain}.int"
                try:
                    os.system(f"external_tools/TMalign {protein_path} {interface_path} -m alignment/matrix.out > alignment/out.tm")
                except:
                    print(f"TM-align did not run for {protein} and {template}_{chain}.")
                parse_tmalign(protein_path, interface_path, protein, template, chain)

    if os.path.exists("alignment/matrix.out"):
        os.system("rm alignment/matrix.out")
    
def parse_tmalign(protein_path, interface_path, protein, template, chain):
    with open("alignment/matrix.out", "r") as matrix_file:
        translation = [0.0, 0.0, 0.0]
        rotationMat = [[0.0, 0.0, 0.0] for _ in range(3)]

        for line in matrix_file:
            tokens = line.strip().split()
            if len(tokens) < 5:
                if line.startswith(" Code"):
                    break
                continue
            try:
                row_index = int(tokens[0])
            except (ValueError, IndexError):
                continue
            # Each matrix row is indicated by 0/1/2 for x/y/z, but original code uses 1/2/3
            # From matrix.out, line numbers are actually 0, 1, 2 in first column. Adjust accordingly.
            # tokens[0]: row (0,1,2); tokens[1]: t[m]; tokens[2-4]: u[m][0:2]
            if row_index in (0, 1, 2):
                translation[row_index] = float(tokens[1])
                rotationMat[row_index][0] = float(tokens[2])
                rotationMat[row_index][1] = float(tokens[3])
                rotationMat[row_index][2] = float(tokens[4])

    with open("alignment/out.tm", "r") as tm_file:
        matchDict = {}
        tmscore1 = 0
        tmscore2 = 0
        for line in tm_file:
            if line.startswith("Aligned length"):
                matchcount = int(line.split("=")[1].split(",")[0])
            elif line.startswith("TM-score"):
                tmscore = float(line.split()[1])
                if line.strip()[-8:-1] == "Chain_1":
                    tmscore1 = tmscore
                if line.strip()[-8:-1] == "Chain_2":
                    tmscore2 = tmscore
            elif line.startswith("(\":\""):
                line = tmOuthnd.readline() #sequence 1
                seq1 = list(line)
                line = tmOuthnd.readline() #matching line
                match = list(line)
                line = tmOuthnd.readline() #sequence 2
                seq2 = list(line)

        seq1ResIDs, seq1chainIDs = extract_chain_and_res_ids("protein", protein_path)
        seq2ResIDs, seq2chainIDs = extract_chain_and_res_ids("interface", interface_path)
        index1 = 0
        index2 = 0    
        for s, i in zip(match, range(len(match))):
            if s == ":" or s == ".": 
                seq1String =  seq1chainIDs[index1] + "." + seq1[i] + "." + seq1ResIDs[index1]
                seq2String =  seq2chainIDs[index2] + "." + seq2[i] + "." + seq2ResIDs[index2] # Reference Model
                matchDict[seq2String] = seq1String
            if seq1[i] != "-":
                index1 += 1
            if seq2[i] != "-":
                index2 += 1

    multi_dict = [matchcount, translation, rotationMat, matchDict, max(tmscore1, tmscore2)]
    pickle.dump(multi_dict, open(f"alignment/{protein}_{template}_{chain}.pkl", "wb"), protocol=2)

def extract_chain_and_res_ids(name, path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(name, path)
    resIDs = []
    chainIDs = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    resIDs.append(str(residue.id[1]))
                    chainIDs.append(chain.id)
    return resIDs, chainIDs