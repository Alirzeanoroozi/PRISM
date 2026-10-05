#!/usr/bin/env python

import os,ConfigParser # pyright: ignore[reportMissingImports]
# import MySQLdb as mdb
import pickle

class StructuralAligner:
    global pdbList,templateList,interfacePath,multiprotOutPath,multiprot
    #constructor for structural aligner requires pdbList,templateList and workPath
    def __init__(self,pdbList,templateList,workPath):
        #if os.path.exists("%s/transformation" % workPath):
        #    return
        # fetching data base info
        #db_f = open("../db_info/dbinfo.inc", "r")
        #db_f.readline()
        #my_host = db_f.readline().split("'")[1]
        #my_user = db_f.readline().split("'")[1]
        #my_pass = db_f.readline().split("'")[1]
        #my_db = db_f.readline().split("'")[1]
        #db_f.close()
        # self.con = mdb.connect(host=my_host, user=my_user, passwd=my_pass, db=my_db)
        # self.cur = self.con.cursor()
        self.pdbList = pdbList
        self.templateList = templateList
        currentPath = os.getcwd()
        os.chdir(workPath)
        config = ConfigParser.ConfigParser() #reads prism.ini file to get configuration datas
        config.read('prism.ini')

        self.interfacePath = config.get('Structural_Alignment','interface_path')
        self.tmalign = config.get('External_Tools','tmalign') #where executable located.
        if not (os.path.exists("alignment")):
            os.mkdir("alignment",0777)
        self.aligner()
        
        if os.path.exists("matrix.out"):
            os.system("rm matrix.out")
        os.chdir(currentPath)

    def aligner (self):
        for protein in self.pdbList:
            self.alignIndividual(protein)
        # if self.con:
        #     self.con.close()

    #def alignIndividual(self,protein):
    #    for interface in self.templateList:
    #        interface = interface[0:6]
    #        check = os.path.exists(self.interfacePath+"/%s_%s.int" % (interface,interface[4])) and os.path.exists(self.interfacePath+"/%s_%s.int" % (interface,interface[5]))
    #        if check:
    #            self.checkTMAlign(protein,interface)
    def alignIndividual(self, protein):
        for interface in self.templateList:
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


    #def checkTMAlign(self,protein,interface):
    #    # check = (protein == "pdb1") or (protein == "pdb2")
    #    # if check:
    #    self.runLocal(protein, interface, interface[4])
    #    self.runLocal(protein, interface, interface[5])
    #    # else:
    #    #     try: 
    #    #         try:
    #    #             self.cur.execute("SELECT content FROM alignment_result where interface=%s && chain=%s && target=%s && method='TM-align'",(interface,interface[4],protein))
    #    #             row = self.cur.fetchone()
    #    #             if row is None:
    #    #                 self.runTMalign(protein,interface,interface[4])
    #    #         except:
    #    #             print "Failed to fetch DB results for alignment_result alignment between %s and %s_%s" % (protein,interface,interface[4])
    #    #             self.runTMalign(protein,interface,interface[4])
    #    #         try:
    #    #             self.cur.execute("SELECT content FROM alignment_result where interface=%s && chain=%s && target=%s && method='TM-align'",(interface,interface[5],protein))
    #    #             row = self.cur.fetchone()
    #    #             if row is None:
    #    #                 self.runTMalign(protein,interface,interface[5])
    #    #         except:
    #    #             print "Failed to fetch DB results for alignment_result alignment between %s and %s_%s" % (protein,interface,interface[5])
    #    #             self.runTMalign(protein,interface,interface[5])
    #    #     except:
    #    #         print "Failed to run TM-align alignment between %s and %s_%s" % (protein,interface,interface)

    #def runTMalign(self, protein, interface, chain):  # this funtion puts tmalign result to mysql table
    #    # name should be changed afterwards
    #    proteinPath = "surfaceExtract/" + protein + ".asa.pdb"
    #    interfaceSide = self.interfacePath + "/%s_%s.int" % (interface, chain)
    #    if os.path.exists(proteinPath) and os.path.exists(interfaceSide):
    #        consoleOut = "alignment/out.tm"
    #        try:
    #            os.system(self.tmalign + " %s %s -m matrix.out > %s" % (proteinPath, interfaceSide, consoleOut))
    #        except:
    #            print
    #            "TM-align did not run for %s and %s." % (interfaceSide, proteinPath)
    #        multiDict = self.parseTMalign(proteinPath, interfaceSide)
    #        # try:
    #        #     self.cur.execute(
    #        #         "SELECT content FROM alignment_result WHERE interface=%s AND chain=%s AND target=%s AND method='TM-align'",
    #        #         (interface, chain, protein))
    #        #     row = self.cur.fetchone()
    #        #     if row is None:
    #        #         self.cur.execute(
    #        #             "INSERT INTO alignment_result (interface,chain,target,method,content) VALUES (%s,%s,%s,%s,%s)",
    #        #             (interface, chain, protein, "TM-align", pickle.dumps(multiDict, protocol=2)))
    #        #     else:
    #        #         self.cur.execute(
    #        #             "UPDATE alignment_result SET content=%s WHERE interface=%s AND chain=%s AND target=%s AND method='TM-align'",
    #        #             (pickle.dumps(multiDict, protocol=2), interface, chain, protein, "TM-align"))
    #        #     self.con.commit()
    #        # except mdb.Error, e:
    #        #     print
    #        #     e
    #        #     self.con.rollback()

    def runLocal(self,protein,interface,chain):#if pdb is provided from user do not put it into mysql
        proteinPath = "surfaceExtract/"+protein + ".asa.pdb"
        interfaceSide = self.interfacePath + "/%s_%s.int" % (interface,chain)
        if os.path.exists(proteinPath) and os.path.exists(interfaceSide):
            consoleOut = "alignment/out.tm"
            try:
                os.system(self.tmalign + " %s %s  -m matrix.out > %s" % (proteinPath,interfaceSide,consoleOut))
            except:
                print "TM-align did not run for %s and %s." % (interfaceSide,proteinPath)
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
                tokens = line.split()
                if len(tokens) < 5:
                    continue
                try:
                    row = int(tokens[0])
                except ValueError:
                    continue
                if row not in (0, 1, 2):
                    continue
                translation[row] = float(tokens[1])
                rotationMat[row][0] = float(tokens[2])
                rotationMat[row][1] = float(tokens[3])
                rotationMat[row][2] = float(tokens[4])

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
                    tmscore = float(line.split()[1].strip())
                    if "Chain_1" in line:
                        tmscore1 = tmscore
                    if "Chain_2" in line:
                        tmscore2 = tmscore
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


