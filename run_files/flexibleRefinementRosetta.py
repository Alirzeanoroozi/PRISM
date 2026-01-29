#!/usr/bin/env python

import os
import configparser
from .fiberdockInterfaceExtractor import fiberdockInterfaceExtractor as fb

class FlexibleRefinement:
    def __init__(self, workPath, jobId):
        try:
            self.workPath = workPath
            self.jobId = jobId
            self.sizes = {}
            self.ca_Num = {}
            self.lastResidue = {}
            self.aDict = {}
            self.currentPath = os.getcwd()
            os.chdir(workPath)
            config = configparser.ConfigParser()  # reads configuration data
            config.read('prism.ini')
            self.rosettaprepack = config.get('External_Tools', 'rosettaprepack')
            self.rosettadock = config.get('External_Tools', 'rosettadock')
            self.rosettadb = config.get('External_Tools', 'rosettadb')
            self.rosettaIntScoreThreshold = config.getfloat('Flexible_Refinement', 'rosettaIntScoreThreshold')
            self.interfacePath = config.get('Structural_Alignment', 'interface_path')
            self.pdbPath = os.path.abspath(config.get('Pdb_Folder', 'pdb_path'))  # defines where to download pdbs wrt workpath and changes into abs
            self.out_folder = "%s/%s" % (config.get('Results_storage', 'output_folder'), self.jobId)
            self.global_ros = "../../%s" % (self.out_folder)

            if not os.path.exists("../../" + config.get('Results_storage', 'output_folder')):
                print ("creating output directory\n")
                os.makedirs("../../" + config.get('Results_storage', 'output_folder'), exist_ok=True)

            if not os.path.exists("flexibleRefinement"):
                os.makedirs("flexibleRefinement", exist_ok=True)
            if not os.path.exists("flexibleRefinement/energies"):
                os.makedirs("flexibleRefinement/energies", exist_ok=True)
            if not os.path.exists("flexibleRefinement/structures"):
                os.makedirs("flexibleRefinement/structures", exist_ok=True)
        except Exception as e:
            print("Exception during initialization: %s" % str(e))
            os.chdir(self.currentPath)

    def refiner(self):
        try:
            energy_Structure = []
            energyFile = "flexibleRefinement/energies/refinement_energies"
            filehndOut = open(energyFile, "w")
            if not os.path.exists("transformation/passedFiles"):
                print("passedFiles does not exists, check other steps.")
            else:
                filehnd = open("transformation/passedFiles", "r")
                for passed in filehnd.readlines():
                    passed = passed.split()
                    passed0 = passed[0].strip()
                    passed1 = passed[1].strip()

                    # totalscore: total energy and should be negative
                    # intscore: Rosetta I_sc score (interaction energy)
                    # structure: the path to the refined structure
                    totalscore, intscore, structure = self.calculateEnergy(passed0, passed1)

                    if intscore != "-":
                        energy_Structure.append(["%s\t%s\t%s" % (passed0, passed1, intscore), structure])
                        filehndOut.writelines("%s\t%s\t%s\n" % (passed0, passed1, intscore))
                filehnd.close()
                filehndOut.close()

            predFile = "%s/refinement_energies" % (self.global_ros)
            if len(energy_Structure) > 0:
                os.system("cp %s %s/refinement_energies" % (energyFile, self.global_ros))
                os.system("chmod 775 %s" % (predFile))

            os.chdir(self.currentPath)
            return energy_Structure
        except Exception as e:
            print("Exception during refiner: %s" % str(e))
            os.chdir(self.currentPath)
            return []

    def calculateEnergy(self, passed0, passed1):
        try:
            target1ChainNameMapping = {}
            target2ChainNamesMapping = {}
            original_chainsT1 = list(self.getChainsFromPdb("transformation/" + passed0))
            original_chainsT2 = list(self.getChainsFromPdb("transformation/" + passed1))
            all_characters = [chr(ord("A") + x) for x in range(26)] + [chr(ord("a") + x) for x in range(26)] + [chr(ord("0") + x) for x in range(10)]
            if len(original_chainsT1) + len(original_chainsT2) > 62:
                print("The number of chains in target1 and target2 exceeds 62, can't assign unique IDs for chains. Rosetta run skipped")
                return "-", "-", "-"

            new_chainsT1 = all_characters[:len(original_chainsT1)]
            new_chainsT2 = all_characters[len(original_chainsT1):len(original_chainsT1)+len(original_chainsT2)]
            chainReNamesMapping = {}
            for i in range(len(original_chainsT1)):
                target1ChainNameMapping[original_chainsT1[i]] = new_chainsT1[i]
                chainReNamesMapping[new_chainsT1[i]] = original_chainsT1[i]
            for i in range(len(original_chainsT2)):
                target2ChainNamesMapping[original_chainsT2[i]] = new_chainsT2[i]
                chainReNamesMapping[new_chainsT2[i]] = original_chainsT2[i]

            self.changeChainName("transformation/" + passed0, target1ChainNameMapping)
            self.changeChainName("transformation/" + passed1, target2ChainNamesMapping)

            combinedpath = self.combinePdb(passed0, passed1)
            partnerChains = "".join(new_chainsT1) + "_" + "".join(new_chainsT2)
            jPath = os.getcwd()

            os.system("%s -database %s -s %s -partners %s -ex1 -ex2aro -out:file:scorefile %s/%s_prepeck_score.sc -overwrite -ignore_zero_occupancy false -detect_disulf false" % (self.rosettaprepack, self.rosettadb, combinedpath, partnerChains, jPath + "/flexibleRefinement/energies/", combinedpath))
            prepackedfile = combinedpath.split(".pdb")[0] + "_0001.pdb"
            os.system("mv %s flexibleRefinement/" % (prepackedfile.split("/")[1]))
            os.system("%s -database %s -s %s -docking_local_refine -partners %s -ex1 -ex2aro -overwrite -ignore_zero_occupancy false -detect_disulf false -out:path:score %s" % (self.rosettadock, self.rosettadb, prepackedfile, partnerChains, jPath + "/flexibleRefinement/energies/"))
            outName = (combinedpath.split("/")[1]).split(".pdb")[0]
            outPdb = outName + "_0001_0001.pdb"

            if not os.path.exists("%s" % (self.global_ros)):
                os.makedirs(self.global_ros, exist_ok=True)
                os.system("chmod 775 %s" % (self.global_ros))

            totalscore = '-'
            interactionscore = '-'

            if os.path.exists("flexibleRefinement/energies/score.sc"):
                os.system("mv flexibleRefinement/energies/score.sc flexibleRefinement/energies/%s_score.sc" % (outName))
                os.system("mv %s flexibleRefinement/structures/%s" % (outPdb, outPdb))

                self.changeChainName("flexibleRefinement/structures/" + outPdb, chainReNamesMapping)
                scorefile = open("flexibleRefinement/energies/%s_score.sc" % (outName), "r")

                for i, line in enumerate(scorefile):
                    if i == 2:
                        temp = line.split()
                        totalscore = float(temp[1].strip())
                        interactionscore = float(temp[5].strip())
                        break
                scorefile.close()
            else:
                print("Could not find a score file for " + outPdb)

            structureList = {0: [], 1: []}
            rosettaOutStructure = "flexibleRefinement/structures/%s" % (outPdb)
            if os.path.exists(rosettaOutStructure) and interactionscore != '-' and interactionscore <= self.rosettaIntScoreThreshold:
                filehnd = open(rosettaOutStructure, "r")
                for chain in original_chainsT1:
                    for l in filehnd:
                        if l[:3] == "TER":
                            break
                        if l[:4] == "ATOM" and l[21] == chain:
                            structureList[0].append(l)
                print("length of struct 0: " + str(len(structureList[0])))

                for chain in original_chainsT2:
                    for l in filehnd:
                        if l[:3] == "TER":
                            break
                        if l[:4] == "ATOM" and l[21] == chain:
                            structureList[1].append(l)
                print("length of struct 1: " + str(len(structureList[1])))

                filehnd.close()

                temp1 = passed0.split("_")
                temp2 = passed1.split("_")
                key1 = temp1[0]+","+temp1[2].split(".")[0]+","+temp2[2].split(".")[0]+","+str(interactionscore)

                check = True
                if key1 in self.aDict:
                    check = False
                else:
                    self.aDict[key1] = 1
                try:
                    os.system("cp %s %s" % (rosettaOutStructure, os.path.join(self.global_ros, outPdb)))
                except Exception as e:
                    print(e)

                intResPath = os.path.join(self.global_ros, "%s.intRes.txt" % (outPdb))
                try:
                    fb.fiberdockInterfaceExtractor(os.path.join(self.global_ros, outPdb), intResPath, structureList[0], structureList[1])
                except Exception as e:
                    print(e)
                    print("Rosetta output visualization and interface residues file were not created.")

                return totalscore, str(interactionscore), os.path.join(self.out_folder, outPdb)
            else:
                if not os.path.exists(rosettaOutStructure):
                    print("structure file couldn't be found %s!!" % (rosettaOutStructure))
                elif interactionscore == '-':
                    print("interactionscore is '-' for structure %s!!" % (rosettaOutStructure))
                elif interactionscore > self.rosettaIntScoreThreshold:
                    print("interactionscore %s exceeds threshold for structure %s!!" % (interactionscore, rosettaOutStructure))
            return "-", "-", "-"
        except Exception as e:
            print("Exception during calculateEnergy: %s" % str(e))
            return "-", "-", "-"

    def getChainsFromPdb(self, pdbPath):
        try:
            chains = set([])
            pdbfile = open(pdbPath, "r")
            for line in pdbfile:
                if line[0:4] == "ATOM":
                    chains.add(line[21])
            pdbfile.close()
            return chains
        except Exception as e:
            print("Exception during getChainsFromPdb: %s" % str(e))
            return set([])

    def changeChainName(self, proteinPath, chainNameMapping):
        try:
            prevfile = open(proteinPath, "r")
            newfile = open(proteinPath + "_temp", "w")
            for line in prevfile:
                if len(line) >= 4 and line[0:4] == "ATOM":
                    if line[21] in chainNameMapping:
                        line = line[0:21] + chainNameMapping[line[21]] + line[22:]
                newfile.write(line)
            prevfile.close()
            newfile.close()
            os.system("cp %s %s" % (proteinPath + "_temp", proteinPath))
        except Exception as e:
            print("Exception during changeChainName: %s" % str(e))

    def combinePdb(self, passed0, passed1):
        try:
            ppath0 = "transformation/" + passed0
            ppath1 = "transformation/" + passed1
            passed0list = passed0.split("_")
            passed1list = passed1.split("_")
            combinedpath = "flexibleRefinement/" + passed0list[0] + "_" + passed0list[2].split(".")[0] \
                           + "_" + passed0list[3] + "_" + passed1list[2].split(".")[0] + "_" + passed1list[3] \
                           + ".rosetta.pdb"
            p0file = open(ppath0, "r")
            p1file = open(ppath1, "r")
            combinedfile = open(combinedpath, "w")
            for line in p0file:
                if line[0:4] == "ATOM":
                    combinedfile.write(line)
            combinedfile.write("TER\n")
            for line in p1file:
                if line[0:4] == "ATOM":
                    combinedfile.write(line)
            combinedfile.write("END")
            p0file.close()
            p1file.close()
            combinedfile.close()
            return combinedpath
        except Exception as e:
            print("Exception during combinePdb: %s" % str(e))
            return ""
