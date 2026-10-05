#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""irmsd_new_fixed.py

I-RMSD implementation based on:
Mendez et al. PROTEINS 52:51-67 (2003)

This is a modernized, Biopython-3.12+/Biopython-1.8+ compatible version of the
original script. The key fix is replacing the removed Bio.SubsMat module with
Bio.Align.substitution_matrices, while keeping Bio.pairwise2 for minimal code
changes.

USAGE
-----
python irmsd_new_fixed.py \
  <model_complex.pdb> <model_receptor_chains> <model_ligand_chains> \
  <ref_complex.pdb>   <ref_receptor_chains>   <ref_ligand_chains>

Examples:
  python irmsd_new_fixed.py model.pdb A B ref.pdb A B
  python irmsd_new_fixed.py model.pdb AB C ref.pdb AB C

Notes:
- Each receptor/ligand chain argument may be a *string* of chain IDs.
  Example: receptor is chains A and B => pass "AB".
- This script expects each complex to be in ONE PDB file that contains BOTH
  receptor and ligand chains. If you have separate unbound receptor/ligand PDBs
  (DB5.5), you must merge them into a single complex-like file per structure.

"""

import sys
import math
import warnings
import itertools
from typing import Dict, List, Tuple

import numpy as np

from Bio.PDB import PDBParser
from Bio import pairwise2
from Bio.Align import substitution_matrices


# -----------------------------
# Substitution matrix handling
# -----------------------------

def load_blosum62_as_dict() -> Dict[Tuple[str, str], float]:
    """Return BLOSUM62 as a dict-of-pairs for pairwise2.align.globalds.

    pairwise2.align.globalds expects a dict with keys like ('A','C') -> score.
    Modern Biopython provides matrices via Bio.Align.substitution_matrices.
    """
    m = substitution_matrices.load("BLOSUM62")
    alphabet = list(m.alphabet)
    md: Dict[Tuple[str, str], float] = {}
    for i, a in enumerate(alphabet):
        for j, b in enumerate(alphabet):
            md[(a, b)] = float(m[i, j])
    return md


BLOSUM62 = load_blosum62_as_dict()


# -----------------------------
# Residue / sequence utilities
# -----------------------------

def getAA_Alphabet(resnam: str) -> str:
    """Map 3-letter residue names to 1-letter codes for standard amino acids."""
    aa = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }
    return aa.get(resnam, '')


def parseChainResiduesFromStructure(modelFile: str, chain: str):
    """Parse a PDB file and return (residue_list, residue_sequence) for a chain.

    Only residues that have backbone atoms CA, N, C, O are kept.
    """
    parser = PDBParser(QUIET=True)
    modelStructure = parser.get_structure('pdbStruct', modelFile)

    resList = []
    resSeq = ""

    # Structure -> Model(0) -> Chain -> Residue
    if chain not in modelStructure[0]:
        raise KeyError(f"Chain '{chain}' not found in {modelFile}")

    for res in modelStructure[0][chain]:
        if res.has_id('CA') and res.has_id('N') and res.has_id('C') and res.has_id('O'):
            resList.append(res)
            resSeq += getAA_Alphabet(res.get_resname())

    return resList, resSeq


# -----------------------------
# Symmetry / chain permutation
# -----------------------------

def getChainIndices(currChainsList: List[str], actualChainString: str) -> List[int]:
    return [actualChainString.index(ch) for ch in currChainsList]


def almostIdentical(seq1: str, seq2: str) -> bool:
    """Treat chains as symmetric if global alignment has <=6 mismatches.

    Uses BLOSUM62, gap open -10, gap extend -1 (original settings).
    """
    matrix = BLOSUM62
    alns = pairwise2.align.globalds(seq1, seq2, matrix, -10, -1)
    if not alns:
        return False

    a0, a1 = alns[0][0], alns[0][1]
    numMismatches = sum(1 for i in range(len(a0)) if a0[i] != a1[i])
    return numMismatches <= 6


def getSymmetricChainsList(modelFile: str, modelChains: str):
    """Return all permutations for groups of symmetric chains.

    If no symmetric chains detected, returns [modelChains] (single string).
    Else returns list of chain-strings like ['AB','BA', ...].
    """
    symmetricChainsList: List[List[str]] = []
    symmetricChainsStrings: List[str] = []

    chainSequences: Dict[str, str] = {}
    chainsDict: Dict[str, List[str]] = {}

    # Step 1: sequences for each chain
    for ch in modelChains:
        _, chSeq = parseChainResiduesFromStructure(modelFile, ch)
        chainSequences[ch] = chSeq
        chainsDict[ch] = [ch]

    # Step 2: group identical chains
    for indx1 in range(len(modelChains)):
        ch1 = modelChains[indx1]
        for indx2 in range(indx1 + 1, len(modelChains)):
            ch2 = modelChains[indx2]
            if ch1 in chainSequences and ch2 in chainSequences:
                if almostIdentical(chainSequences[ch1], chainSequences[ch2]):
                    chainsDict[ch1].extend(chainsDict[ch2])
                    del chainsDict[ch2]

    # Step 3: count permutations
    numSymmChainStrings = 0
    for ch in chainsDict:
        if len(chainsDict[ch]) > 1:
            numSymmChainStrings += math.factorial(len(chainsDict[ch]))

    if numSymmChainStrings == 0:
        return [modelChains]

    for _ in range(numSymmChainStrings):
        symmetricChainsList.append(list(itertools.repeat('0', len(modelChains))))

    # Step 4: fill permutations
    for ch in chainsDict:
        currSymmList = chainsDict[ch]
        currIndicesList = getChainIndices(currSymmList, modelChains)

        sclcount = 0
        while sclcount < numSymmChainStrings:
            for perm in itertools.permutations(currSymmList):
                for i in range(len(currIndicesList)):
                    indx = currIndicesList[i]
                    symmetricChainsList[sclcount][indx] = perm[i]
                sclcount += 1

    # Step 5: list -> strings
    for chlist in symmetricChainsList:
        symmetricChainsStrings.append(''.join(chlist))

    return symmetricChainsStrings


# -----------------------------
# Alignment-based residue pruning
# -----------------------------

def getCommonResiduesList(refFile: str, refChains: str, modelFile: str, modelChains: str):
    """Align corresponding chains and remove residues missing in either structure.

    Returns:
        (refComponentResidues, modelComponentResidues)
    where each is a list of Residue objects in the *same order* and *same count*.
    """
    refComponent = []
    modelComponent = []

    matrix = BLOSUM62

    for refch, mdlch in zip(refChains, modelChains):
        modelResidueList, modelSequence = parseChainResiduesFromStructure(modelFile, mdlch)
        refResidueList, refSequence = parseChainResiduesFromStructure(refFile, refch)

        alns = pairwise2.align.globalds(modelSequence, refSequence, matrix, -10, -1)
        if not alns:
            continue

        modelAligned = alns[0][0]
        refAligned = alns[0][1]

        modelResiduesToRemove = []
        refResiduesToRemove = []

        mcount = 0
        rcount = 0

        for i in range(len(modelAligned)):
            if modelAligned[i] != '-' and refAligned[i] != '-':
                rcount += 1
                mcount += 1
            elif modelAligned[i] == '-':
                # model gap => remove ref residue
                if rcount < len(refResidueList):
                    refResiduesToRemove.append(refResidueList[rcount])
                rcount += 1
            elif refAligned[i] == '-':
                # ref gap => remove model residue
                if mcount < len(modelResidueList):
                    modelResiduesToRemove.append(modelResidueList[mcount])
                mcount += 1

        for itm in modelResiduesToRemove:
            if itm in modelResidueList:
                modelResidueList.remove(itm)
        for itm in refResiduesToRemove:
            if itm in refResidueList:
                refResidueList.remove(itm)

        refComponent.extend(refResidueList)
        modelComponent.extend(modelResidueList)

    return refComponent, modelComponent


# -----------------------------
# Interface detection
# -----------------------------

def combineCoords(receptorInterface: Dict[int, np.ndarray], ligandInterface: Dict[int, np.ndarray]) -> np.ndarray:
    coords: List[np.ndarray] = []
    coords = [receptorInterface[indx] for indx in sorted(receptorInterface.keys())]
    coords.extend([ligandInterface[indx] for indx in sorted(ligandInterface.keys())])
    return np.array(coords)


def getInterfaceResidues(refReceptorResList, refLigandResList, modelReceptorResList, modelLigandResList):
    """Return coordinates (CA) of interface residues based on reference interface."""

    interDis = 10.0
    safeDis = 20.0

    refRecInterface: Dict[int, np.ndarray] = {}
    refLigInterface: Dict[int, np.ndarray] = {}

    modelRecInterface: Dict[int, np.ndarray] = {}
    modelLigInterface: Dict[int, np.ndarray] = {}

    for rrescount, rres in enumerate(refReceptorResList):
        for lrescount, lres in enumerate(refLigandResList):
            ca_dist = lres['CA'] - rres['CA']

            if ca_dist > safeDis:
                continue
            elif ca_dist < interDis:
                if rrescount not in refRecInterface:
                    refRecInterface[rrescount] = rres['CA'].get_coord()
                    modelRecInterface[rrescount] = modelReceptorResList[rrescount]['CA'].get_coord()
                if lrescount not in refLigInterface:
                    refLigInterface[lrescount] = lres['CA'].get_coord()
                    modelLigInterface[lrescount] = modelLigandResList[lrescount]['CA'].get_coord()
                continue
            else:
                gotonext = False
                for ratm in rres:
                    if gotonext:
                        break
                    for latm in lres:
                        if (ratm - latm) < interDis:
                            if rrescount not in refRecInterface:
                                refRecInterface[rrescount] = rres['CA'].get_coord()
                                modelRecInterface[rrescount] = modelReceptorResList[rrescount]['CA'].get_coord()
                            if lrescount not in refLigInterface:
                                refLigInterface[lrescount] = lres['CA'].get_coord()
                                modelLigInterface[lrescount] = modelLigandResList[lrescount]['CA'].get_coord()
                            gotonext = True
                            break

    referenceInterface = combineCoords(refRecInterface, refLigInterface)
    modelInterface = combineCoords(modelRecInterface, modelLigInterface)

    return referenceInterface, modelInterface


# -----------------------------
# RMSD (Kabsch-like / original)
# -----------------------------

def calcRotMatAndEigenValues(a: np.ndarray, b: np.ndarray, mu: np.ndarray):
    U = np.zeros((3, 3))
    for i in range(3):
        U = U + np.outer(b[:, i], a[:, i])

    detU = np.linalg.det(U)

    # Reflection check
    if np.abs(detU + 1.0) < 0.001:
        indx = int(mu.argmin(0))
        mu[indx] = -mu[indx]
        U = U - 2 * np.outer(b[:, indx], a[:, indx])

    return U, mu


def getRMSD(ra: np.ndarray, rb: np.ndarray) -> float:
    """Compute RMSD between two Nx3 coordinate arrays."""

    if ra.shape != rb.shape:
        raise ValueError(f"Coordinate arrays must have same shape. Got {ra.shape} vs {rb.shape}")

    lenSuperpose = ra.shape[0]
    if lenSuperpose == 0:
        return 0.0

    # Translation
    centroidA = ra.sum(axis=0) / lenSuperpose
    centroidB = rb.sum(axis=0) / lenSuperpose

    ra = ra - centroidA
    rb = rb - centroidB

    # R matrix
    R = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            R[i, j] = (rb[:, i] * ra[:, j]).sum()

    Rt = R.conj().transpose()
    eigvals, eigvecs = np.linalg.eig(np.dot(Rt, R))

    # Numerical safety
    eigvals = np.where(eigvals < 0, 0, eigvals)

    mu = np.sqrt(eigvals)
    a = eigvecs

    b = np.zeros((3, 3))
    for k in range(3):
        if mu[k] == 0:
            b[:, k] = 0
        else:
            b[:, k] = (np.dot(R, a[:, k])) / mu[k]

    U, mu = calcRotMatAndEigenValues(a, b, mu)

    Dsquare = np.power(ra, 2).sum() + np.power(rb, 2).sum() - 2.0 * (mu.sum())

    if Dsquare > 0.0:
        rmsd = math.sqrt(Dsquare) / math.sqrt(lenSuperpose)
    else:
        rmsd = 0.0

    return float(rmsd)


# -----------------------------
# Main IRMSD routine
# -----------------------------

def getIRMSD(argv: List[str]) -> float:
    """Compute least IRMSD over symmetric chain permutations."""

    # Silence BiopythonDeprecationWarning from pairwise2
    warnings.filterwarnings('ignore')

    if len(argv) != 7:
        raise SystemExit(
            "\n".join([
                "ERROR: wrong number of arguments.",
                "Usage:",
                f"  {argv[0]} <model.pdb> <modelR> <modelL> <ref.pdb> <refR> <refL>",
                "Example:",
                f"  {argv[0]} model.pdb A B ref.pdb A B",
            ])
        )

    modelFile = argv[1]
    modelRch = argv[2]
    modelLch = argv[3]

    refFile = argv[4]
    refRch = argv[5]
    refLch = argv[6]

    # Symmetric chain combinations
    if len(modelRch) > 1:
        mdlReceptorChainCombinations = getSymmetricChainsList(modelFile, modelRch)
    else:
        mdlReceptorChainCombinations = [modelRch]

    if len(modelLch) > 1:
        mdlLigandChainCombinations = getSymmetricChainsList(modelFile, modelLch)
    else:
        mdlLigandChainCombinations = [modelLch]

    least_irmsd = 100000.0

    for mdlRchain in mdlReceptorChainCombinations:
        for mdlLchain in mdlLigandChainCombinations:
            refReceptor, modelReceptor = getCommonResiduesList(refFile, refRch, modelFile, mdlRchain)
            refLigand, modelLigand = getCommonResiduesList(refFile, refLch, modelFile, mdlLchain)

            refCoords, modelCoords = getInterfaceResidues(refReceptor, refLigand, modelReceptor, modelLigand)

            if len(refCoords) == 0 or len(modelCoords) == 0:
                curr_irmsd = 100.0
            else:
                curr_irmsd = getRMSD(modelCoords, refCoords)

            if curr_irmsd < least_irmsd:
                least_irmsd = curr_irmsd

    return least_irmsd


def main():
    val = getIRMSD(sys.argv)
    print(f"{val:.3f}")


if __name__ == '__main__':
    main()