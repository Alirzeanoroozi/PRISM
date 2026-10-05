import os,ConfigParser
#a class to do preprocessing on the proteins before surface extraction
#Written by Alper Baspinar
class PreProcessor:
    global leftTarget,rightTarget,pdbPath,currentPath,amino #global PDBlist and pdbPath to reach from everywhere in the class
    #constructor requires pdbList and a workPath to run
    def __init__(self,leftTarget,rightTarget,workPath):

        # ---- DEFINE ALL ATTRIBUTES FIRST (CRITICAL FIX) ----
        self.leftTarget = leftTarget
        self.rightTarget = rightTarget
        self.currentPath = os.getcwd()
        self.workPath = workPath

        self.pdbPath = "pdb"
        self.amino = [
            'ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU',
            'MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR'
        ]

        # ---- ORIGINAL LOGIC BELOW (UNCHANGED) ----
        os.chdir(workPath)

        config = ConfigParser.ConfigParser()
        config.read('prism.ini')

        if not os.path.exists("preprocess"):
            os.mkdir("preprocess", 0777)

        if os.path.exists("%s/surfaceExtract" % workPath):
            pass#return
        
    def prepareProtein(self):    
        leftList = []
        rightList = []
        if len(self.leftTarget) == len(self.rightTarget):
            for index in range(len(self.leftTarget)):
                proteinLeft = self.leftTarget[index]
                proteinRight = self.rightTarget[index]
                path1 = self.pdbPath+"/%s.pdb" % proteinLeft[0:4]
                path2 = self.pdbPath+"/%s.pdb" % proteinRight[0:4]
                check1 = False
                check2 = False 
                if os.path.exists(path1): #checks if protein exists in the pdb path or not
                    extra = proteinLeft[4:] #if pdb also contains chains pdb files should be splitted
                    if len(extra) == 0:
                        residueList = self.splitNotChain(path1)
                    else:
                        residueList = self.splitChains(extra,path1)
                    if self.chainWriter(proteinLeft,residueList):
                        check1 = True
                else:
                    continue
                
                if os.path.exists(path2): #checks if protein exists in the pdb path or not
                    extra = proteinRight[4:] #if pdb also contains chains pdb files should be splitted
                    if len(extra) == 0:
                        residueList = self.splitNotChain(path2)
                    else:
                        residueList = self.splitChains(extra,path2)
                    if self.chainWriter(proteinRight,residueList):
                        check2 = True
                else:
                    continue
            
                if check1 and check2:
                    leftList.append(proteinLeft)
                    rightList.append(proteinRight)             
        os.chdir(self.currentPath) #changes running path to previous one    
        return leftList,rightList

    def splitNotChain(self,path):
        residueList = []
        alternateLocationDic = {}
        icodeDic = {}
        filehnd = open(path,"r")
        for line in filehnd.readlines():
            if line[:3] == "END":
                break
            if line[0:4] == "ATOM":
                chain = line[21]
                atomName = line[12:16]
                resName = line[17:20]
                residueSeq = line[22:26]
                altLoc = line[16] #icode and alternate location handled
                icode = line[26]
                key = chain+residueSeq
                key2 = key+atomName
                check1 = alternateLocationDic.has_key(key2)
                check2 = icodeDic.has_key(key)
                check3 = resName in self.amino
                #these lines are needed to detect alternateLocation or icode differences, the first one is taken...
                if check1 == False and check2 == False and check3 == True:
                    alternateLocationDic[key2] = altLoc
                    icodeDic[key] = icode
                    residueList.append(line)
                elif check1 == False and check2 == True and check3 == True:
                    alternateLocationDic[key2] = altLoc
                    if icodeDic[key] == icode:
                        residueList.append(line)
        filehnd.close()
        return residueList
        
        
    def splitChains(self,extra,path):
        residueList = []
        alternateLocationDic = {}
        icodeDic = {}
        filehnd = open(path,"r")
        for line in filehnd.readlines():
            if line[:3] == "END":
                break
            if line[0:4] == "ATOM":
                chain = line[21]
                if chain in extra:
                    atomName = line[12:16]
                    residueSeq = line[22:26]
                    resName = line[17:20]
                    altLoc = line[16]
                    icode = line[26]
                    key = chain+residueSeq
                    key2 = key+atomName
                    check1 = alternateLocationDic.has_key(key2)
                    check2 = icodeDic.has_key(key)
                    check3 = resName in self.amino
                    #these lines are needed to detect alternateLocation or icode differences, the first one is taken...
                    if check1 == False and check2 == False and check3 == True:
                        alternateLocationDic[key2] = altLoc
                        icodeDic[key] = icode
                        residueList.append(line)
                    elif check1 == False and check2 == True and check3 == True:
                        alternateLocationDic[key2] = altLoc
                        if icodeDic[key] == icode:
                            residueList.append(line)
        filehnd.close()
        return residueList
        
    #writes to a file
    def chainWriter(self,protein,residueList):
        writePath = "preprocess/"+protein+".pdb"
        if len(residueList) != 0:
            newFile = open(writePath,"w")
            for residue in residueList:
                newFile.writelines(residue)
            newFile.writelines("END")
            newFile.close()
            return True
        else:
            return False

