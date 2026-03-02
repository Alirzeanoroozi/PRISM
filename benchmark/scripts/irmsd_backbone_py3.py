#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
I-rmsd as defined in:
Mendez et al.
PROTEINS: Structure, Function and Genetics 52:51-67 (2003)
Assessment of Blind Predictions of Protein-Protein Interactions:
Current Status of Docking Methods.

Backbone-atom version adapted for Python 3.
"""

import argparse
import math
import itertools
import warnings

import numpy
from Bio.PDB import PDBParser
from Bio.Align import PairwiseAligner, substitution_matrices

from utills import getAA_Alphabet


def parseChainResiduesFromStructure(modelFile, chain):
    parser = PDBParser(QUIET=True)
    modelStructure = parser.get_structure("pdbStruct", modelFile)
    resList = []
    resSeq = ""
    for res in modelStructure[0][chain]:
        if res.has_id("CA") and res.has_id("N") and res.has_id("C") and res.has_id("O"):
            resList.append(res)
            resSeq = resSeq + getAA_Alphabet(res.get_resname())
    return resList, resSeq


def _make_aligner():
    matrix = substitution_matrices.load("BLOSUM62")
    aligner = PairwiseAligner()
    aligner.substitution_matrix = matrix
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -1
    return aligner


def getSymmetricChainsList(modelFile, modelChains):
    symmetricChainsList = []
    symmetricChainsStrings = []

    chainSequences = {}
    chainsDict = {}

    for ch in modelChains:
        _, chSeq = parseChainResiduesFromStructure(modelFile, ch)
        chainSequences[ch] = chSeq
        chainsDict[ch] = [ch]

    for indx1 in range(len(modelChains)):
        ch1 = modelChains[indx1]
        for indx2 in range(indx1 + 1, len(modelChains)):
            ch2 = modelChains[indx2]
            if ch1 in chainSequences and ch2 in chainSequences:
                if almostIdentical(chainSequences[ch1], chainSequences[ch2]):
                    chainsDict[ch1].extend(chainsDict[ch2])
                    del chainsDict[ch2]

    numSymmChainStrings = 0
    for ch in chainsDict:
        if len(chainsDict[ch]) > 1:
            numSymmChainStrings = numSymmChainStrings + math.factorial(len(chainsDict[ch]))

    if numSymmChainStrings == 0:
        return [modelChains]

    for _ in range(numSymmChainStrings):
        symmetricChainsList.append(list(itertools.repeat(0, len(modelChains))))

    for ch in chainsDict:
        currSymmList = chainsDict[ch]
        currIndicesList = getChainIndices(currSymmList, modelChains)

        sclcount = 0
        while sclcount < numSymmChainStrings:
            for perm in itertools.permutations(currSymmList):
                for i in range(len(currIndicesList)):
                    indx = currIndicesList[i]
                    symmetricChainsList[sclcount][indx] = perm[i]
                sclcount = sclcount + 1

    for chlist in symmetricChainsList:
        chlist = "".join(chlist)
        symmetricChainsStrings.append(chlist)

    return symmetricChainsStrings


def getChainIndices(currChainsList, actualChainString):
    return [actualChainString.index(ch) for ch in currChainsList]


def almostIdentical(seq1, seq2):
    aligner = _make_aligner()
    alignments = aligner.align(seq1, seq2)
    alignment = alignments[0]
    aln0, aln1 = alignment[0], alignment[1]
    numMismatches = sum(1 for i in range(len(aln0)) if aln0[i] != aln1[i])
    return numMismatches <= 6


def getCommonResiduesList(refFile, refChains, modelFile, modelChains):
    refComponent = []
    modelComponent = []
    aligner = _make_aligner()

    for refch, mdlch in zip(refChains, modelChains):
        modelResidueList, modelSequence = parseChainResiduesFromStructure(modelFile, mdlch)
        refResidueList, refSequence = parseChainResiduesFromStructure(refFile, refch)
        alignments = aligner.align(modelSequence, refSequence)
        alignment = alignments[0]
        modelAligned = alignment[0]
        refAligned = alignment[1]

        modelResiduesToRemove = []
        refResiduesToRemove = []
        mcount = 0
        rcount = 0

        for i in range(len(modelAligned)):
            if modelAligned[i] != "-" and refAligned[i] != "-":
                rcount += 1
                mcount += 1
            elif modelAligned[i] == "-":
                refResiduesToRemove.append(refResidueList[rcount])
                rcount += 1
            elif refAligned[i] == "-":
                modelResiduesToRemove.append(modelResidueList[mcount])
                mcount += 1

        for itm in modelResiduesToRemove:
            modelResidueList.remove(itm)
        for itm in refResiduesToRemove:
            refResidueList.remove(itm)

        refComponent.extend(refResidueList)
        modelComponent.extend(modelResidueList)

    return refComponent, modelComponent


def getBackboneCoord(res):
    return [res["N"].get_coord(), res["CA"].get_coord(), res["C"].get_coord(), res["O"].get_coord()]


def getInterfaceResidues(refReceptorResList, refLigandResList, modelReceptorResList, modelLigandResList):
    """Returns a list of interface residues for a structure (backbone atoms)."""
    interDis = 10.0
    safeDis = 20.0
    refRecInterface = {}
    refLigInterface = {}
    modelRecInterface = {}
    modelLigInterface = {}

    for rrescount in range(len(refReceptorResList)):
        rres = refReceptorResList[rrescount]
        for lrescount in range(len(refLigandResList)):
            lres = refLigandResList[lrescount]
            if lres["CA"] - rres["CA"] > safeDis:
                continue
            if lres["CA"] - rres["CA"] < interDis:
                if rrescount not in refRecInterface:
                    refRecInterface[rrescount] = getBackboneCoord(rres)
                    modelRecInterface[rrescount] = getBackboneCoord(modelReceptorResList[rrescount])
                if lrescount not in refLigInterface:
                    refLigInterface[lrescount] = getBackboneCoord(lres)
                    modelLigInterface[lrescount] = getBackboneCoord(modelLigandResList[lrescount])
                continue

            gotonext = False
            for ratm in rres:
                if gotonext:
                    break
                for latm in lres:
                    if ratm - latm < interDis:
                        if rrescount not in refRecInterface:
                            refRecInterface[rrescount] = getBackboneCoord(rres)
                            modelRecInterface[rrescount] = getBackboneCoord(modelReceptorResList[rrescount])
                        if lrescount not in refLigInterface:
                            refLigInterface[lrescount] = getBackboneCoord(lres)
                            modelLigInterface[lrescount] = getBackboneCoord(modelLigandResList[lrescount])
                        gotonext = True
                        break

    referenceInterface = combineBackboneCoords(refRecInterface, refLigInterface)
    modelInterface = combineBackboneCoords(modelRecInterface, modelLigInterface)
    return referenceInterface, modelInterface


def combineBackboneCoords(receptorInterface, ligandInterface):
    coords = []
    for indx in sorted(receptorInterface.keys()):
        coords.extend(receptorInterface[indx])
    for indx in sorted(ligandInterface.keys()):
        coords.extend(ligandInterface[indx])
    return numpy.array(coords)


def getRMSD(ra, rb):
    """ra and rb are Nx3 arrays."""
    R = numpy.zeros((3, 3), dtype=float)
    b = numpy.zeros((3, 3), dtype=float)
    lenSuperpose = ra.shape[0]
    centroidA = ra.sum(axis=0) / lenSuperpose
    centroidB = rb.sum(axis=0) / lenSuperpose
    ra = ra - centroidA
    rb = rb - centroidB
    for i in range(3):
        for j in range(3):
            R[i, j] = (rb[:, i] * ra[:, j]).sum()
    Rtranspose = R.transpose()
    eig_vals, a = numpy.linalg.eigh(numpy.dot(Rtranspose, R))
    eig_vals = numpy.maximum(eig_vals, 0.0)
    mu = numpy.sqrt(eig_vals)
    for k in range(3):
        if mu[k] == 0.0:
            continue
        b[:, k] = (numpy.dot(R, a[:, k])) / mu[k]
    _, mu = calcRotMatAndEigenValues(a, b, mu)
    Dsquare = numpy.power(ra, 2).sum() + numpy.power(rb, 2).sum() - 2.0 * (mu.sum())
    if Dsquare > 0.0:
        rmsd = math.sqrt(Dsquare) / math.sqrt(lenSuperpose)
    else:
        rmsd = 0.0
    return rmsd


def calcRotMatAndEigenValues(a, b, mu):
    U = numpy.zeros((3, 3), dtype=float)

    for i in range(3):
        U = U + numpy.outer(b[:, i], a[:, i])

    detU = numpy.linalg.det(U)
    if numpy.abs(detU + 1.0) < 0.001:
        _, indx = mu.min(0), mu.argmin(0)
        mu[indx] = -mu[indx]
        U = U - 2 * numpy.outer(b[:, indx], a[:, indx])

    return U, mu


def _parse_args():
    parser = argparse.ArgumentParser(
        description="Compute interface RMSD using all backbone atoms.",
    )
    parser.add_argument("model_file", help="Model PDB file")
    parser.add_argument("model_receptor_chains", help="Model receptor chain string, e.g. AB")
    parser.add_argument("model_ligand_chains", help="Model ligand chain string, e.g. C")
    parser.add_argument("reference_file", help="Reference PDB file")
    parser.add_argument("reference_receptor_chains", help="Reference receptor chain string")
    parser.add_argument("reference_ligand_chains", help="Reference ligand chain string")
    return parser.parse_args()


def getIRMSD():
    warnings.filterwarnings("ignore")
    args = _parse_args()
    modelFile = args.model_file
    modelRch = args.model_receptor_chains
    modelLch = args.model_ligand_chains
    refFile = args.reference_file
    refRch = args.reference_receptor_chains
    refLch = args.reference_ligand_chains

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
    print(f"{least_irmsd:.3f}")


if __name__ == "__main__":
    print("iRMSD with all backbone atoms")
    getIRMSD()
