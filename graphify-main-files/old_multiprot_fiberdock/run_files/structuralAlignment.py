#!/usr/bin/env python
#Written by Alper Baspinar

import os,ConfigParser
#CLI/fiberdock port: MySQL multiprot cache removed; all alignments written to files
import pickle

class StructuralAligner:
    global pdbList,templateList,interfacePath,multiprotOutPath,multiprot
    #constructor for structural aligner requires pdbList,templateList and workPath
    def __init__(self,pdbList,templateList,workPath):
    	if os.path.exists("%s/transformation" % workPath):
    		return
        self.pdbList = pdbList
        self.templateList = templateList
        currentPath = os.getcwd()
        os.chdir(workPath)
        config = ConfigParser.ConfigParser() #reads prism.ini file to get configuration datas
        config.read('prism.ini')
        self.interfacePath = config.get('Structural_Alignment','interface_path')
        self.multiprot = config.get('External_Tools','multiprot') #where executable located.
        
        if not (os.path.exists("alignment")):
		os.mkdir("alignment",0777)
        self.aligner()
        if os.path.exists("2_sets.res"):
        	os.system("rm 2_sets.res")
        if os.path.exists("log_multiprot.txt"):
        	os.system("rm log_multiprot.txt")
        os.chdir(currentPath)

    def aligner (self):
        for protein in self.pdbList:
                self.alignIndividual(protein)

    def alignIndividual(self,protein):
        for interface in self.templateList:
            interface = interface[0:6]
            check = os.path.exists(self.interfacePath+"/%s_%s.int" % (interface,interface[4])) and os.path.exists(self.interfacePath+"/%s_%s.int" % (interface,interface[5]))
	    if check:
            	self.checkMultiProt(protein,interface)

    def checkMultiProt(self,protein,interface):
        #CLI/fiberdock port: no DB cache -- always run multiprot and write result to file
        self.runLocal(protein,interface,interface[4])
        self.runLocal(protein,interface,interface[5])

    def runLocal(self,protein,interface,chain):#run multiprot and store result as a pickle file
	proteinPath = "surfaceExtract/"+protein + ".asa.pdb"
        interfaceSide = self.interfacePath + "/%s_%s.int" % (interface,chain)
        if os.path.exists(proteinPath) and os.path.exists(interfaceSide):
            consoleOut = "alignment/%s.multiprot" % (protein)
            try:
                os.system(self.multiprot + " %s %s > %s" % (interfaceSide,proteinPath,consoleOut))
            except:
                print "Multiprot did not run for %s and %s." % (interfaceSide,proteinPath)
            multiDict = self.parseMultiprot()
   	    fileName = "alignment/%s_%s_%s" % (interface,chain,protein)
	    filehnd = open(fileName,"w")
	    filehnd.write(pickle.dumps(multiDict,protocol=2))
	    filehnd.close() 	
		

    def parseMultiprot(self): #parses three results of the multiprot
	if os.path.exists("2_sol.res"):
		multiOuthnd = open("2_sol.res",'r')
		multiDict = {}
		count = 0
		line = multiOuthnd.readline()
		while not(line == ""):
		    if count == 3:
			break
		    elif line[:12] == "Solution Num":
			multiOuthnd.readline() #empty string
			line = multiOuthnd.readline() #Mult Corres Score : ?
			line = line.split(":")
			matchcount = int(line[1])
			line = multiOuthnd.readline() #Reference Molecule : ?
			line = line.split(":")
			refMol = int(line[1])
			multiOuthnd.readline() # Molecule
			line = multiOuthnd.readline() # Trans ? ? ? ? ? ?
			line = (line.split(":")[1]).strip().split()
			transV = [float(line[0]),float(line[1]),float(line[2]),float(line[3]),float(line[4]),float(line[5])]
			multiOuthnd.readline() # RMSD
			multiOuthnd.readline() # empty string
			multiOuthnd.readline() # Match List
			matchDict = {}
			for m in range(matchcount):
			    line = multiOuthnd.readline()
			    if line[:3] == "End":
				break
			    else:
				line = line.strip()
				line = line.split()
				matchDict[line[refMol]] = line[1-refMol] #refMol can be 0 or 1 we always want refMol = 0 which is interface

			multiDict[count] = [matchcount,refMol,transV,matchDict] # multiDict contains all information in the multiprot file no need to read multiprot file again
			count += 1
			line = multiOuthnd.readline()
		    else:
			line = multiOuthnd.readline()
		multiOuthnd.close()
		os.system("rm 2_sol.res")
		return multiDict
	else:
		return -1

