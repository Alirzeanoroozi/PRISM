#!/usr/bin/env python
#Written by Alper Baspinar
import string,os,ConfigParser  

#uses pops to extract surface of the proteins
class SurfaceExtractor:
    global RSATHRESHOLD, SCFFTHRESHOLD, pdbList, external_tool_choice, external_tool,currentPath
    #constructor of the class requires pdbList and workPath to run
    def __init__(self, pdbList, workPath):

        # ---- DEFINE ALL ATTRIBUTES FIRST (CRITICAL FIX) ----
        self.pdbList = pdbList
        self.currentPath = os.getcwd()
        self.workPath = os.path.abspath(workPath)
        self.preprocessPath = os.path.join(self.workPath, "preprocess")
        self.surfacePath = os.path.join(self.workPath, "surfaceExtract")



        self.RSATHRESHOLD = None
        self.SCFFTHRESHOLD = None
        self.external_tool_choice = None
        self.external_tool = None

        # ---- ORIGINAL LOGIC BELOW (UNCHANGED) ----
        #if os.path.exists("%s/alignment" % workPath):
        #    return
        self.skip_surface = os.path.exists(os.path.join(workPath, "alignment"))


        #os.chdir(workPath)
        os.chdir(self.workPath)


        config = ConfigParser.ConfigParser()
        config.read('prism.ini')

        self.RSATHRESHOLD = config.getfloat('Surface_Extractor', 'rsathreshold')
        self.SCFFTHRESHOLD = config.getfloat('Surface_Extractor', 'scffthreshold')
        self.external_tool_choice = config.getint('Surface_Extractor', 'external_tool_choice')
        if self.external_tool_choice == 0:
            self.external_tool = config.get('External_Tools', 'pops')
        elif self.external_tool_choice == 1:
            self.external_tool = config.get('External_Tools', 'naccess')
        else:
            os.chdir(self.currentPath)
            return

        #if not os.path.exists("surfaceExtract"):
        #    os.mkdir("surfaceExtract", 0777)
        if not os.path.exists(self.surfacePath):
            os.makedirs(self.surfacePath, 0777)
#################################### START SURFACE EXTRACTION METHODS #####################################
    def surfaceExtractor(self):
        pdbList = []
        for protein in self.pdbList:
            if self.findSurface(protein):
                pdbList.append(protein)
        os.chdir(self.currentPath) 
        return pdbList

    def findSurface(self,protein):
        success = self.runExternal(protein)
        if success == 1:
            asaDict = self.extractSurface(protein)
            return self.outputAsaFile(protein,asaDict)
        else:
            return False
 
    def runExternal(self, proteinname):
        proteinPath = "preprocess/"+proteinname + ".pdb" #protein path
        if self.external_tool_choice == 0:
            popsOutputFile = "surfaceExtract/"+proteinname + ".rsa" #protein result path
            errorFile = "surfaceExtract/error"

            if os.path.exists(proteinPath):
                os.system("%s --pdb %s --residueOut --popsOut %s > %s" % (self.external_tool,proteinPath, popsOutputFile, errorFile))
            else:
                print "protein %s does not exist in the protein path, check previous stages " % proteinname
                return -1

            if os.path.exists(errorFile) == False:
                print "pops did not run for %s" % proteinname
                return -1

            errorhandler = open(errorFile)
            errorstring = errorhandler.read()
            if string.find(errorstring,"STOP") != -1 or string.find(errorstring,"IEEE") != -1 or string.find(errorstring,"Error") != -1 or string.find(errorstring,"Clean termination") == -1: #checks whether pops successfully ended or not
                print "pops error! " + proteinname
                return -1
            errorhandler.close()
            if os.path.exists("sigma.out"):
                os.system("rm sigma.out")
            return 1
        else:
            naccessOut = "surfaceExtract/naccessOut"
            naccessError = "surfaceExtract/naccessError"

            if os.path.exists(proteinPath):
                os.system("%s %s > %s 2>%s" % (self.external_tool,proteinPath,naccessOut,naccessError))
            else:
                print "protein %s does not exist in the protein path, check previous stages " % proteinname
                return -1

            if os.path.exists(naccessError) == False:
                print "naccess did not run for %s" % proteinname
                return -1

            errorhandler = open(naccessError)
            errorstring = errorhandler.read()
            if string.find(errorstring,"STOP") != -1 or string.find(errorstring,"IEEE") != -1 or string.find(errorstring,"error") != -1:
                print "naccess error! " + proteinname
                return -1
            errorhandler.close()
            rsafile = proteinname + ".rsa"
            asafile = proteinname + ".asa"
            logfile = proteinname + ".log"
            naccessOutputFile = "surfaceExtract/"+proteinname + ".rsa" #protein result path
            if os.path.exists(rsafile):
                os.system("mv %s %s"%(rsafile,naccessOutputFile))
            if os.path.exists(asafile):
                os.system("rm %s" % (asafile))
            if os.path.exists(logfile):
                os.system("rm %s" % (logfile))
            return 1

    def extractSurface(self,proteinname): #reads pops data together with pdb file to get asa 
        rsaPath = "surfaceExtract/"+proteinname + ".rsa"
        proteinPath = "preprocess/%s.pdb" % (proteinname) 
        rsahandler = open(rsaPath,"r")
        asahandler = open(proteinPath,"r")
        rsalist = []
        rsaline = ""
        asaDict = {}
        if self.external_tool_choice == 0:
            for rsaline in rsahandler.readlines():
                standard = self.StandardData(rsaline[:3])
                if standard != -1:
                    rsatemp = rsaline.strip().split()
                    resseq = rsatemp[2]
                    reschn = rsatemp[1]
                    try:
                        absoluteacc = float(rsatemp[5])
                    except:
                        print "pops output handle error for protein %s" % proteinname
                        return asaDict
                    relativeacc = absoluteacc*100/standard

                    if relativeacc > self.RSATHRESHOLD:
                        rsalist.append((reschn, resseq))
            rsahandler.close()
        else:
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
                            print "naccess output handle error for protein %s" % proteinname
                            return asaDict
                        relativeacc = absoluteacc*100/standard
                        if relativeacc > self.RSATHRESHOLD:
                            rsalist.append((reschn,resseq))
            rsahandler.close()

        for asaline in asahandler.readlines():
            if asaline[:3] == "END":
                break
            if asaline[:4] == "ATOM":
                atom = asaline[13:15]
                resNo = asaline[22:26].strip() 
                chain = asaline[21]
                serialNumber = asaline[6:11]
                for reschn,resseq in rsalist:
                    if atom == "CA" and resNo == resseq and chain == reschn:
                        asaDict[serialNumber] = asaline

        asahandler.close()

        return self.findScaffold(asaDict, proteinPath)

    def findScaffold(self, asaDict, proteinPath):
        pdbhnd = open(proteinPath, "r")
        #find nearby residues and add to the asalist
        for line in pdbhnd.readlines():
            if line[:3] == "END":
                break
            serialNumber = line[6:11]
            if line[:4] == "ATOM" and line[13:15] == "CA" and (serialNumber in asaDict.keys()) == False:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                for key in asaDict.keys():
                    xcoor = float(asaDict[key][30:38])
                    ycoor = float(asaDict[key][38:46])
                    zcoor = float(asaDict[key][46:54])
                    dist = ((xcoor-x)**2 + (ycoor-y)**2 + (zcoor-z)**2)**0.5
                    if dist <= self.SCFFTHRESHOLD:
                        if asaDict.has_key(serialNumber) == False:
                            asaDict[serialNumber] = line

        pdbhnd.close()
        return asaDict


    def StandardData(self,resname):
        residueDict = {'ALA':107.95,'CYS':134.28,'ASP':140.39,'GLU':172.25,'PHE':199.48,'GLY':80.1,'HIS':182.88,'ILE':175.12,\
                       'LYS':200.81,'LEU':178.63,'MET':194.15,'ASN':143.94,'PRO':136.13,'GLN':178.5,'ARG':238.76,'SER':116.5,\
                       'THR':139.27,'VAL':151.44,'TRP':249.36,'TYR':212.76}
        
        if residueDict.has_key(resname):
            return residueDict[resname]               
        else: 
            return -1
    
    
    #writes the output file protein.asa.pdb 
    def outputAsaFile(self,protein, asaDict):
        if len(asaDict) != 0:
            asahnd = open("surfaceExtract/%s.asa.pdb" % (protein), "w")
            for key in sorted(asaDict.keys()):
                asahnd.writelines(asaDict[key])
            asahnd.writelines("END")
            asahnd.close()
            return True
        else:
            return False
        
