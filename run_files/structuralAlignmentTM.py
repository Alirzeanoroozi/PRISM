import os
import configparser
import pickle

class StructuralAligner:
    def __init__(self, queries, templates):

        self.queries = queries
        self.templates = templates
        config = configparser.ConfigParser() # reads prism.ini file to get configuration datas
        config.read('prism.ini')

        self.interfacePath = config.get('Structural_Alignment','interface_path')
        self.tmalign = config.get('External_Tools','tmalign') # where executable located.
        
        os.makedirs("alignment", exist_ok=True)
        for protein in self.queries:
            self.alignIndividual(protein)
        
        if os.path.exists("matrix.out"):
            os.system("rm matrix.out")

    def alignIndividual(self, protein):
        for interface in self.templates:
            interface = interface[0:6]

            chainA = interface[4]
            chainB = interface[5]

            alnA = os.path.join("alignment", "%s_%s_%s" % (interface, chainA, protein))
            alnB = os.path.join("alignment", "%s_%s_%s" % (interface, chainB, protein))

            intA = self.interfacePath + "/%s_%s.int" % (interface, chainA)
            intB = self.interfacePath + "/%s_%s.int" % (interface, chainB)

            # Skip interface completely if interface files are missing

            if not (os.path.exists(intA) and os.path.exists(intB)):
                #print('Skipping alignment', interface, 'due to missing files.')
                continue

            # Run only missing alignments
            if not os.path.exists(alnA):
                self.runLocal(protein, interface, chainA)

            if not os.path.exists(alnB):
                self.runLocal(protein, interface, chainB)

    def runLocal(self, protein, interface, chain):#if pdb is provided from user do not put it into mysql
        proteinPath = "surfaceExtract/" + protein + ".asa.pdb"
        interfaceSide = self.interfacePath + "/%s_%s.int" % (interface,chain)
        if os.path.exists(proteinPath) and os.path.exists(interfaceSide):
            consoleOut = "alignment/out.tm"
            try:
                os.system(self.tmalign + " %s %s  -m matrix.out > %s" % (proteinPath,interfaceSide,consoleOut))
            except:
                print("TM-align did not run for %s and %s." % (interfaceSide,proteinPath))
            multiDict = self.parseTMalign(proteinPath, interfaceSide)
            fileName = "alignment/%s_%s_%s" % (interface,chain,protein)
            filehnd = open(fileName,"w")
            filehnd.write(pickle.dumps(multiDict,protocol=2))
            filehnd.close()     
        
    def parseTMalign(self, proteinPath, interfaceSide): #parses three results of the TM align file
        if os.path.exists("matrix.out") and os.path.exists("alignment/out.tm"):
            matrixHnd = open("matrix.out",'r')
            tmOuthnd = open("alignment/out.tm",'r')
            multiDict = {}
            matchcount = 0
            refMol = 0 # our reference is always the interface side (complex B in TM files)
            translation = [0, 0, 0]
            rotationMat = [[0 for x in range(3)] for y in range(3)] 
            for line in matrixHnd:
                if line[1] == "1":
                    line = line.split()
                    translation[0] = float(line[1])
                    rotationMat[0][0] = float(line[2])
                    rotationMat[0][1] = float(line[3])
                    rotationMat[0][2] = float(line[4])
                elif line[1] == "2":
                    line = line.split()
                    translation[1] = float(line[1])
                    rotationMat[1][0] = float(line[2])
                    rotationMat[1][1] = float(line[3])
                    rotationMat[1][2] = float(line[4])
                elif line[1] == "3":
                    line = line.split()
                    translation[2] = float(line[1])
                    rotationMat[2][0] = float(line[2])
                    rotationMat[2][1] = float(line[3])
                    rotationMat[2][2] = float(line[4])
                elif line[:5] == " Code":
                    break

            transV = [translation, rotationMat]
            
            matchDict = {}
            tmscore1 = 0
            tmscore2 = 0
            line = tmOuthnd.readline()
            while not(line == ""):
                if line[:14] == "Aligned length":
                    line = line.split("=")
                    line = line[1].split(",")
                    matchcount = int(line[0])
                elif line[:8] == "TM-score": # we take into account the TMScore with higher value
                    if line.strip()[-8:-1] == "Chain_1":
                        tmscore1 = line.split()[1].strip()
                        tmscore1 = float(tmscore1)
                    if line.strip()[-8:-1] == "Chain_2":
                        tmscore2 = line.split()[1].strip()
                        tmscore2 = float(tmscore2)
                elif line[:4] == "(\":\"":
                    line = tmOuthnd.readline() #sequence 1
                    seq1 = list(line)
                    line = tmOuthnd.readline() #matching line
                    match = list(line)
                    line = tmOuthnd.readline() #sequence 2
                    seq2 = list(line)

                    # get res IDs for seq1 and seq2
                    pdb1hnd = open(proteinPath, 'r')
                    pdb2hnd = open(interfaceSide, 'r')
                    seq1ResIDs = []
                    seq2ResIDs = []
                    seq1chainIDs = []
                    seq2chainIDs = []
                    for line1 in pdb1hnd:
                        if line1[:3] == "END":
                            break
                        if(line1[:4] == "ATOM" and line1[13:15] == "CA"):
                            seq1ResIDs.append(line1[22:26].strip())
                            seq1chainIDs.append(line1[21])
                    for line2 in pdb2hnd.readlines():
                        if line2[:3] == "END":
                            break
                        if(line2[:4] == "ATOM" and line2[13:15] == "CA"):
                            seq2ResIDs.append(line2[22:26].strip())
                            seq2chainIDs.append(line2[21])
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
                line = tmOuthnd.readline()
            matrixHnd.close()
            tmOuthnd.close()
            if matchcount == 0 :
                return -1 
            multiDict[0] = [matchcount,refMol,transV,matchDict,max(tmscore1,tmscore2)] # multiDict contains all information in the multiprot file no need to read multiprot file again
            return multiDict
        else:
            return -1


