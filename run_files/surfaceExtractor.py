# Written by Alper Baspinar
# uses pops to extract surface of the proteins

import os
import configparser

class SurfaceExtractor:
    def __init__(self, queries):
        self.queries = queries

        config = configparser.ConfigParser() #reads prism.ini file for configuration datas
        config.read('prism.ini')
        self.RSATHRESHOLD = config.getfloat('Surface_Extractor','rsathreshold')
        self.SCFFTHRESHOLD = config.getfloat('Surface_Extractor','scffthreshold')
        self.external_tool = config.get('External_Tools','naccess')

    def surfaceExtractor(self):
        for protein in self.queries:
            if self.runExternal(protein) == 1:
                asaDict = self.extractSurface(protein)
                self.outputAsaFile(protein, asaDict)

    def runExternal(self, proteinname):
        protein_path = "preprocess/" + proteinname + ".pdb" #protein path
        os.makedirs("surfaceExtract", exist_ok=True)
        naccessOut = "surfaceExtract/naccessOut"
        naccessError = "surfaceExtract/naccessError"

        if os.path.exists(protein_path):
            os.system("%s %s > %s 2>%s" % (self.external_tool, protein_path, naccessOut, naccessError))
        else:
            print("protein %s does not exist in the protein path, check previous stages " % proteinname)
            return -1

        if os.path.exists(naccessError) == False:
            print("naccess did not run for %s" % proteinname)
            return -1

        errorhandler = open(naccessError)
        errorstring = errorhandler.read()
        if "STOP" in errorstring or "IEEE" in errorstring or "Error" in errorstring:
            print("naccess error! " + proteinname)
            return -1
        errorhandler.close()
        
        rsa_file = proteinname + ".rsa"
        asa_file = proteinname + ".asa"
        logfile = proteinname + ".log"
        naccess_output_file = "surfaceExtract/" + proteinname + ".rsa" #protein result path
        if os.path.exists(rsa_file):
            os.system("mv %s %s"%(rsa_file, naccess_output_file))
        if os.path.exists(asa_file):
            os.system("rm %s" % (asa_file))
        if os.path.exists(logfile):
            os.system("rm %s" % (logfile))
        return 1

    def extractSurface(self, proteinname): #reads pops data together with pdb file to get asa 
        rsaPath = "surfaceExtract/" + proteinname + ".rsa"
        proteinPath = "preprocess/%s.pdb" % (proteinname) 
        rsalist = []
        asaDict = {}

        with open(rsaPath, "r") as rsahandler, open(proteinPath, "r") as asahandler:
            for rsaline in rsahandler.readlines():
                if rsaline[:3] == "RES":
                    rsaline = rsaline.strip()
                    resR = rsaline[4:7]
                    standard = self.StandardData(resR)
                    if standard != -1:
                        reschn = rsaline[8]
                        resseq = rsaline[9:13].strip()
                        try:
                            absoluteacc = float(rsaline[14:22])
                        except:
                            print("naccess output handle error for protein %s" % proteinname)
                            return asaDict
                        relativeacc = absoluteacc * 100 / standard
                        if relativeacc > self.RSATHRESHOLD:
                            rsalist.append((reschn,resseq))

            for asaline in asahandler.readlines():
                if asaline[:3] == "END":
                    break
                if asaline[:4] == "ATOM":
                    atom = asaline[13:15]
                    resNo = asaline[22:26].strip() 
                    chain = asaline[21]
                    serialNumber = asaline[6:11]
                    for reschn, resseq in rsalist:
                        if atom == "CA" and resNo == resseq and chain == reschn:
                            asaDict[serialNumber] = asaline

        return self.findScaffold(asaDict, proteinPath)

    def findScaffold(self, asaDict, proteinPath):
        with open(proteinPath, "r") as pdbhnd:
            for line in pdbhnd.readlines():
                if line[:3] == "END":
                    break
                serialNumber = line[6:11]
                if line[:4] == "ATOM" and line[13:15] == "CA" and serialNumber not in asaDict:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    for key in asaDict:
                        xcoor = float(asaDict[key][30:38])
                        ycoor = float(asaDict[key][38:46])
                        zcoor = float(asaDict[key][46:54])
                        dist = ((xcoor-x)**2 + (ycoor-y)**2 + (zcoor-z)**2)**0.5
                        if dist <= self.SCFFTHRESHOLD:
                            if serialNumber not in asaDict:
                                asaDict[serialNumber] = line

        return asaDict

    def StandardData(self, resname):
        residueDict = {'ALA':107.95,'CYS':134.28,'ASP':140.39,'GLU':172.25,'PHE':199.48,'GLY':80.1,'HIS':182.88,'ILE':175.12,\
                       'LYS':200.81,'LEU':178.63,'MET':194.15,'ASN':143.94,'PRO':136.13,'GLN':178.5,'ARG':238.76,'SER':116.5,\
                       'THR':139.27,'VAL':151.44,'TRP':249.36,'TYR':212.76}
        
        return residueDict.get(resname, -1)
    
    # writes the output file protein.asa.pdb 
    def outputAsaFile(self, protein, asaDict):
        if len(asaDict) != 0:
            asahnd = open("surfaceExtract/%s.asa.pdb" % (protein), "w")
            for key in sorted(asaDict.keys()):
                asahnd.writelines(asaDict[key])
            asahnd.writelines("END")
            asahnd.close()
            return True
        else:
            return False
