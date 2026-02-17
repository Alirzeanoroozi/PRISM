#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# This program uses the Biopython distribution. (http://www.biopython.org)

# 11/30/2011
# @Author Shruthi Viswanath
# modified by Attila Gursoy to include all backbone atoms

''' 
I-rmsd as defined here: 
Mendez et al.
PROTEINS: Structure, Function and Genetics 52:51-67 (2003)
Assessment of Blind Predictions of Protein-Protein Interactions: Current Status of Docking Methods.
'''

import os
import sys
import math
import numpy
import itertools
import warnings
from Bio.PDB import PDBParser
from Bio.Align import PairwiseAligner, substitution_matrices

from utills import getAA_Alphabet
    
def parseChainResiduesFromStructure(modelFile, chain):
    parser = PDBParser()
    modelStructure = parser.get_structure('pdbStruct', modelFile)
    resList = []
    resSeq = ""
    for res in modelStructure[0][chain]:
        if res.has_id('CA') and res.has_id('N') and res.has_id('C') and res.has_id('O'):
            resList.append(res)
            resSeq = resSeq + getAA_Alphabet(res.get_resname())
    return resList, resSeq
    
def getSymmetricChainsList(modelFile, modelChains):
    symmetricChainsList=[] # e.g. [['A','B'],['B','A']] if A and B are identical chains
    symmetricChainsStrings=[] # e.g. [AB,BA]
    
    chainSequences={}
    chainsDict={}
        
    # Step 1. Get the sequences for all chains given and initialize dictionaries of symmetric chains
    for ch in modelChains:
        (chResList,chSeq)=parseChainResiduesFromStructure(modelFile,ch)
        
        chainSequences[ch]=chSeq 

        chainsDict[ch]=[ch]
           
    # Step 2: Get identical chains in a dictionary. e.g. if chains A and B are identical in the string ABC, it returns {'A':['A','B'],'C':['C']}
    for indx1 in range(len(modelChains)):
        ch1=modelChains[indx1]
        
        for indx2 in range(indx1+1,len(modelChains)):
            ch2=modelChains[indx2]
        
        if ch1 in chainSequences and ch2 in chainSequences: # keys exist
        
            if almostIdentical(chainSequences[ch1],chainSequences[ch2]): # sequences identical
                chainsDict[ch1].extend(chainsDict[ch2])
                del chainsDict[ch2]                 
               
    # Step 3: Initialize the strings of chains list
    numSymmChainStrings=0
    for ch in chainsDict:
         if len(chainsDict[ch])>1:
          numSymmChainStrings=numSymmChainStrings+math.factorial(len(chainsDict[ch]))
    
    if numSymmChainStrings==0: # no symmetric chains, return the chain list as it is 
         return([modelChains])
    
    else:
         for count in range(numSymmChainStrings):
            symmetricChainsList.append(list(itertools.repeat(0,len(modelChains))))
    
    # Step 4: Get all the permutations of the new chain lists 
    for ch in chainsDict:
         currSymmList=chainsDict[ch]
    
         currIndicesList=getChainIndices(currSymmList,modelChains)
         
         sclcount=0
         while sclcount<numSymmChainStrings: # loop over each string in the chains list. 
         
            for perm in itertools.permutations(currSymmList): # loop over each permutation of the current group, eg.(A,B) is one permutation, (B,A) is another
             for i in range(len(currIndicesList)): # loop over the length of the string 
                  indx=currIndicesList[i]
                  symmetricChainsList[sclcount][indx]=perm[i]
    
             sclcount=sclcount+1
    
    # Step 5: Convert the list to strings
    for chlist in symmetricChainsList:
           chlist = ''.join(chlist)
           symmetricChainsStrings.append(chlist)        
    
    return(symmetricChainsStrings)

def getChainIndices(currChainsList, actualChainString):
    return [actualChainString.index(ch) for ch in currChainsList]
        
def almostIdentical(seq1, seq2):
    matrix = substitution_matrices.load("BLOSUM62")
    aligner = PairwiseAligner()
    aligner.substitution_matrix = matrix
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -1
    alignments = aligner.align(seq1, seq2)
    alignment = alignments[0]
    aln0, aln1 = alignment[0], alignment[1]
    numMismatches = sum(1 for i in range(len(aln0)) if aln0[i] != aln1[i])
    return numMismatches <= 6
         
def getCommonResiduesList(refFile, refChains, modelFile, modelChains):
    refComponent = []
    modelComponent = []
    matrix = substitution_matrices.load("BLOSUM62")
    aligner = PairwiseAligner()
    aligner.substitution_matrix = matrix
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -1

    for refch, mdlch in zip(refChains, modelChains):
        (modelResidueList, modelSequence) = parseChainResiduesFromStructure(modelFile, mdlch)
        (refResidueList, refSequence) = parseChainResiduesFromStructure(refFile, refch)
        alignments = aligner.align(modelSequence, refSequence)
        alignment = alignments[0]
        modelAligned = alignment[0]
        refAligned = alignment[1]

        modelResiduesToRemove = []
        refResiduesToRemove = []
        mcount = 0
        rcount = 0

        for i in range(len(modelAligned)):
            if modelAligned[i] != '-' and refAligned[i] != '-':
                rcount += 1
                mcount += 1
            elif modelAligned[i] == '-':
                refResiduesToRemove.append(refResidueList[rcount])
                rcount += 1
            elif refAligned[i] == '-':
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
    return [res['N'].get_coord(), res['CA'].get_coord(), res['C'].get_coord(), res['O'].get_coord()]
 
def getInterfaceResidues(refReceptorResList, refLigandResList, modelReceptorResList, modelLigandResList):
    """Returns a list of the interface residues of a given structure (backbone atoms)."""
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
            if lres['CA'] - rres['CA'] > safeDis:
                continue
            elif lres['CA'] - rres['CA'] < interDis:
                if rrescount not in refRecInterface:
                    refRecInterface[rrescount] = getBackboneCoord(rres)
                    modelRecInterface[rrescount] = getBackboneCoord(modelReceptorResList[rrescount])
                if lrescount not in refLigInterface:
                    refLigInterface[lrescount] = getBackboneCoord(lres)
                    modelLigInterface[lrescount] = getBackboneCoord(modelLigandResList[lrescount])
                continue
            else:
                gotonext = False
                for ratm in rres:
                    if not gotonext:
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
                    else:
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

    coords=numpy.array(coords)
    return(coords)         
    
def combineCoords(receptorInterface, ligandInterface):
    coords = []
    coords = [receptorInterface[indx] for indx in sorted(receptorInterface.keys())]
    coords.extend([ligandInterface[indx] for indx in sorted(ligandInterface.keys())])
    
    coords=numpy.array(coords)
                     
    return(coords)         

def getRMSD(ra, rb):
    """ra and rb are Nx3 arrays."""
    R = numpy.zeros((3, 3))
    b = numpy.zeros((3, 3))
    lenSuperpose = ra.shape[0]
    centroidA = ra.sum(axis=0) / lenSuperpose
    centroidB = rb.sum(axis=0) / lenSuperpose
    ra = ra - centroidA
    rb = rb - centroidB
    for i in range(3):
        for j in range(3):
            R[i, j] = (rb[:, i] * ra[:, j]).sum()
    Rtranspose = R.conj().transpose()
    musquare = numpy.linalg.eig(numpy.dot(Rtranspose, R))[0]
    a = numpy.linalg.eig(numpy.dot(Rtranspose, R))[1]
    mu = numpy.sqrt(musquare)
    for k in range(3):
        b[:, k] = (numpy.dot(R, a[:, k])) / mu[k]
    U, mu = calcRotMatAndEigenValues(a, b, mu)
    Dsquare = numpy.power(ra, 2).sum() + numpy.power(rb, 2).sum() - 2.0 * (mu.sum())
    if Dsquare > 0.0:
        rmsd = math.sqrt(Dsquare) / math.sqrt(lenSuperpose)
    else:
        rmsd = 0.00
    return rmsd

def calcRotMatAndEigenValues(a,b,mu):
    U=numpy.zeros((3,3))
        
    for i in range(3):
          U=U+numpy.outer(b[:,i],a[:,i]) 
         
    # Get the determinant of U
    detU=numpy.linalg.det(U)

    # check the determinant of U for reflection. 
    if numpy.abs(detU+1.0000)<0.001:  # Inversion detected 
           minval,indx = mu.min(0),mu.argmin(0)
           mu[indx]=-mu[indx]
           
           U=U-2*numpy.outer(b[:,indx],a[:,indx]) 
        
    return(U,mu)

    
def getIRMSD():
    warnings.filterwarnings('ignore')
    modelFile = sys.argv[1]
    modelRch = sys.argv[2]
    modelLch = sys.argv[3]
    refFile = sys.argv[4]
    refRch = sys.argv[5]
    refLch = sys.argv[6]

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
    print("%.3f" % least_irmsd)


# MAIN
if __name__ == '__main__':
    print("iRMSD with all backbone atoms")
    getIRMSD() 
