# written by Alper Baspinar
import os,math
import configparser
import pickle

class TransformFilter:
    def __init__(self, queries, templates):
        self.queries = queries
        self.templates = templates
        self.templateSize = {}
        
        config = configparser.ConfigParser() #reads configuration datas
        config.read('prism.ini')
        self.MINIMUM_RESIDUE_MATCH_COUNT = config.getfloat('Transformation_Filtering','minimum_residue_match_count')
        self.MINIMUM_RESIDUE_MATCH_PERCENTAGE = config.getfloat('Transformation_Filtering','minimum_residue_match_percentage')
        self.MINIMUM_HOTSPOT_MATCH_NUMBER = config.getfloat('Transformation_Filtering','minimum_hotspot_match_number')
        self.DIFF_PERCENTAGE = config.getfloat('Transformation_Filtering','diff_percentage')
        self.contactPath = config.get('Transformation_Filtering','contactpath')
        self.hotspotPath = config.get('Transformation_Filtering','hotspotdatapath')
        self.multiprotCount = config.getfloat('Transformation_Filtering','multiprotcount')
        self.hotspotCriterion = config.getfloat('Transformation_Filtering','hotspotcriterion')
        self.hotspotCount = config.getfloat('Transformation_Filtering','hotspotCount')
        self.multiprotOutPath = config.get('Structural_Alignment','multiprot_output')
        # Resolve interfacePath robustly: prefer absolute, then workPath-relative, then project-root template/interfaces
        cfg_interface = config.get('Structural_Alignment','interface_path')
        if os.path.isabs(cfg_interface):
            resolved_interface = cfg_interface
        else:
            # prefer path relative to the job workPath
            candidate_wp = os.path.abspath(os.path.join(self.workPath, cfg_interface))
            if os.path.exists(candidate_wp):
                resolved_interface = candidate_wp
            else:
                # try relative to original project root (currentPath)
                candidate_root = os.path.abspath(os.path.join(self.currentPath, cfg_interface))
                if os.path.exists(candidate_root):
                    resolved_interface = candidate_root
                else:
                    # fallback to canonical template/interfaces under project root
                    resolved_interface = os.path.abspath(os.path.join(self.currentPath, 'template', 'interfaces'))
        self.interfacePath = resolved_interface
        print("TransformationFilter: using interfacePath=%s" % (self.interfacePath))
        # Resolve contactPath and hotspotPath relative to workPath if they are not absolute
        cfg_contact = self.contactPath
        if os.path.isabs(cfg_contact):
            resolved_contact = cfg_contact
        else:
            candidate_wp = os.path.abspath(os.path.join(self.workPath, cfg_contact))
            if os.path.exists(candidate_wp):
                resolved_contact = candidate_wp
            else:
                candidate_root = os.path.abspath(os.path.join(self.currentPath, cfg_contact))
                resolved_contact = candidate_root if os.path.exists(candidate_root) else candidate_wp
        self.contactPath = resolved_contact
        print("TransformationFilter: using contactPath=%s" % (self.contactPath))
        cfg_hot = self.hotspotPath
        if os.path.isabs(cfg_hot):
            resolved_hot = cfg_hot
        else:
            candidate_wp = os.path.abspath(os.path.join(self.workPath, cfg_hot))
            if os.path.exists(candidate_wp):
                resolved_hot = candidate_wp
            else:
                candidate_root = os.path.abspath(os.path.join(self.currentPath, cfg_hot))
                resolved_hot = candidate_root if os.path.exists(candidate_root) else candidate_wp
        self.hotspotPath = resolved_hot
        print("TransformationFilter: using hotspotPath=%s" % (self.hotspotPath))
        self.template_Residue_Count = config.getfloat('Transformation_Filtering','template_residue_count')
        self.contact_Count = config.getfloat('Transformation_Filtering','contact_count')
        self.clashing_distance = config.getfloat('Transformation_Filtering','clashing_distance')
        self.max_clashing_count = config.getfloat('Transformation_Filtering','max_clashing_count')
        self.TM_SCORE_THRESHOLD = config.getfloat('Structural_Alignment','tm_score_threshold')

        os.makedirs("transformation", exist_ok=True)

    def transformer(self):
        # machine-readable list (Left\tRight)
        passed_files_path = os.path.join("transformation", "passedFiles")
        readable_path = os.path.join("transformation", "passed_pairs_readable.txt")
        with open(passed_files_path, "w") as filehnd, open(readable_path, "w") as readable_fh:
            readable_fh.write("Readable list of passed transformation pairs (Left -> Right)\n")
            passedInterfaces = []
        for interface in self.templates:
            interface = interface[:6]
            temp = self.filtering(interface)
            passed_pairs, passed_ints = temp[0], temp[1]
            #print("Found %d passed pair(s) and %d passed interface(s) for %s" % (len(passed_pairs), len(passed_ints), interface))
            for passed in passed_pairs:
                left = passed[0]
                right = passed[1]
                line = left+"\t"+right+"\n"
                filehnd.writelines(line)
                readable_fh.write("Pair: Left=%s  ->  Right=%s\n" % (left, right))
                #print("Passed pair: Left=%s  ->  Right=%s" % (left, right))
            for passInt in passed_ints:
                passedInterfaces.append(passInt)
                line = passInt+"\n"
                readable_fh.write("Passed interface: %s\n" % passInt)
                #filehnd2.writelines(line)
        return passedInterfaces
    
    def contactDict(self,interface):
        try:
            contactFile = open(self.contactPath + "/"+interface+'.txt','r')
        except:
            print("Interface %s does not exist..." % interface)
            return -1 #contact file for the interface could not be found!!

        contactDic = {}

        for line in contactFile.readlines():
            contact = line.strip().split()
            key = contact[0]+contact[1]
            contactDic[key] = 1

        contactFile.close()
        return contactDic

    def filtering(self,interface):
        sizePathLeft =  self.interfacePath + "/%s_%s.int" % (interface,interface[4])
        sizePathRight = self.interfacePath + "/%s_%s.int" % (interface,interface[5])

        dContact = self.contactDict(interface)
        if dContact == -1 or (not os.path.exists(sizePathLeft)) or (not os.path.exists(sizePathRight)):
            print('Returning..')
            return [],[]
        else:
            
            #first read interface_chain.int files to calculate size of the template chain
            filehnd = open(sizePathLeft,"r")
            leftKey = "%s_%s" % (interface,interface[4])
            self.templateSize[leftKey] = len(filehnd.readlines())
            filehnd.close()
            filehnd = open(sizePathRight,"r")
            rightKey = "%s_%s" % (interface,interface[5])
            self.templateSize[rightKey] = len(filehnd.readlines())
            #print("Template sizes: %s -> %d, %s -> %d" % (leftKey, self.templateSize[leftKey], rightKey, self.templateSize[rightKey]))

            leftPartner = []
            rightPartner = []
            for index in range(len(self.leftTarget)):
                left1 = "%s_%s_%s" % (interface,interface[4],self.leftTarget[index])
                alignDict1 = self.alignmentDict(self.leftTarget[index],interface,interface[4])
                left2 = "%s_%s_%s" % (interface,interface[4],self.rightTarget[index])
                alignDict2 = self.alignmentDict(self.rightTarget[index],interface,interface[4])
                right1 = "%s_%s_%s" % (interface,interface[5],self.rightTarget[index])
                alignDict3 = self.alignmentDict(self.rightTarget[index],interface,interface[5])
                right2 = "%s_%s_%s" % (interface,interface[5],self.leftTarget[index])
                alignDict4 = self.alignmentDict(self.leftTarget[index],interface,interface[5]) 

                if alignDict1 == -1:
                    print("No alignment for %s (left1)" % left1)
                if alignDict2 == -1:
                    print("No alignment for %s (left2)" % left2)
                if alignDict3 == -1:
                    print("No alignment for %s (right1)" % right1)
                if alignDict4 == -1:
                    print("No alignment for %s (right2)" % right2)

                if alignDict1 != -1 and alignDict3 != -1:
                    #print("Adding partner pair: %s and %s" % (left1, right1))
                    leftPartner.append([left1, alignDict1])
                    rightPartner.append([right1, alignDict3])
                if alignDict2 != -1 and alignDict4 != -1:
                    #print("Adding partner pair: %s and %s" % (left2, right2))
                    leftPartner.append([left2, alignDict2])
                    rightPartner.append([right2, alignDict4])
            leftMatchDict = []
            rightMatchDict = []
            passedInterfaces = []
            for index in range(len(leftPartner)):
                temp1 = self.matchThresholdCheck(leftPartner[index])
                temp2 = self.matchThresholdCheck(rightPartner[index])
                if temp1 != -1:
                    passedInterfaces.append(leftPartner[index][0])
                if temp2 != -1:
                    passedInterfaces.append(rightPartner[index][0])
                if temp1 != -1 and temp2 != -1:
                    leftMatchDict.append([temp1,leftPartner[index][0]])
                    rightMatchDict.append([temp2,rightPartner[index][0]])
    
            passedList = []
            for index in range(len(leftMatchDict)):
                #print('leftMatchDict[index]', leftMatchDict[index][0])
                solutionListLeft = leftMatchDict[index][0][0]
                alignDictLeft = leftMatchDict[index][0][1]
                left = leftMatchDict[index][1]
                solutionListRight = rightMatchDict[index][0][0]
                alignDictRight = rightMatchDict[index][0][1]
                right = rightMatchDict[index][1]
                for solLeft in solutionListLeft:
                      for solRight in solutionListRight:
                        if self.interfaceMatchCheck(dContact,alignDictLeft[solLeft][3],alignDictRight[solRight][3]):
                            leftTarget = self.rotateTarget(left,solLeft,alignDictLeft[solLeft])
                            rightTarget = self.rotateTarget(right,solRight,alignDictRight[solRight])
                            if self.overlap(leftTarget,rightTarget):
                                passedList.append([leftTarget,rightTarget])
            return passedList,passedInterfaces
        
    
    def overlap(self,leftTarget,rightTarget):
        clash_count = 0
        leftTarget = "transformation/"+leftTarget
        rightTarget = "transformation/"+rightTarget
        if not(os.path.exists(leftTarget)):
            print("Left Target %s does not exist." % leftTarget)
            return False
        elif not(os.path.exists(rightTarget)):
            print("Right Target %s does not exist." % rightTarget)
            return False
        else:
            leftCACoor = self.readTarget(leftTarget)
            rightCACoor = self.readTarget(rightTarget)
            for leftcoor in leftCACoor:
                for rightcoor in rightCACoor:
                    if self.getDistance(leftcoor,rightcoor) < self.clashing_distance:
                        clash_count += 1
                        if self.max_clashing_count <= clash_count:
                            return False
            return True

    def readTarget(self,path):
        filehnd = open(path,"r")
        CACoor = []
        for line in filehnd.readlines():
            if line[:3] == "END":
                break
            elif line[:4] == "ATOM" and line[12:16].strip() == "CA":
                x_coor = float(line[30:38])
                y_coor = float(line[38:46])
                z_coor = float(line[46:54])
                coord = [x_coor,y_coor,z_coor]
                CACoor.append(coord)
        filehnd.close()
        return CACoor

    def getDistance(self,coord1,coord2):
        x1,y1,z1 = coord1
        x2,y2,z2 = coord2
        return ((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)**0.5

    def rotateTarget(self,partner,sol,multiList):
        refMol = multiList[1]
        transV = multiList[2]
        temp = partner.split("_")
        template = temp[0]+"_"+temp[1]
        target = temp[2]+".pdb"
        transformPath = template+"_"+target+"_"+str(sol)+"_trans.pdb"
        self.pdbTransform(target,transformPath,transV,refMol)
        return transformPath
######################################################################TRANSFORM SCRIPT#################################
    def pdbTransform(self,target,transformPath,transV,refMol):
        TMalign = False
        if len(transV) == 2:
            x_translate, y_translate, z_translate  = transV[0]
            rotationMatrix = transV[1]
            TMalign = True
        else :
            phi, theta, psi, x_translate, y_translate, z_translate = transV
        
        if refMol == 0:
            if not TMalign:
                rotationMatrix = self.rotationDictExtractor(phi + math.pi, theta + math.pi, psi + math.pi)
            pdbFile = open("preprocess/"+target, 'r')
            transformedPDBFile = open("transformation/"+transformPath, 'w')
            for pdbLine in pdbFile:
                if pdbLine[0:4] == 'ATOM':
                    x = float(pdbLine[30:38].strip())
                    y = float(pdbLine[38:46].strip())
                    z = float(pdbLine[46:54].strip())
                    if not(TMalign):
                        new_x = x*rotationMatrix[0][0] + y*rotationMatrix[1][0] + z*rotationMatrix[2][0] + x_translate
                        new_y = x*rotationMatrix[0][1] + y*rotationMatrix[1][1] + z*rotationMatrix[2][1] + y_translate
                        new_z = x*rotationMatrix[0][2] + y*rotationMatrix[1][2] + z*rotationMatrix[2][2] + z_translate
                    else:
                        new_x = x*rotationMatrix[0][0] + y*rotationMatrix[0][1] + z*rotationMatrix[0][2] + x_translate
                        new_y = x*rotationMatrix[1][0] + y*rotationMatrix[1][1] + z*rotationMatrix[1][2] + y_translate
                        new_z = x*rotationMatrix[2][0] + y*rotationMatrix[2][1] + z*rotationMatrix[2][2] + z_translate
                    pdbLine = '%s%8.3f%8.3f%8.3f%s' %(pdbLine[0:30], new_x, new_y, new_z, pdbLine[54:len(pdbLine)])
                transformedPDBFile.write(pdbLine)
            pdbFile.close()
            transformedPDBFile.close()
        elif refMol == 1:
            if not TMalign:
                rotationMatrix = self.transposeRotationDictExtractor(phi + math.pi, theta + math.pi, psi + math.pi)
            pdbFile = open("preprocess/"+target, 'r')
            transformedPDBFile = open("transformation/"+transformPath, 'w')
            for pdbLine in pdbFile:
                if pdbLine[0:4] == 'ATOM':
                    x = float(pdbLine[30:38].strip()) - x_translate
                    y = float(pdbLine[38:46].strip()) - y_translate
                    z = float(pdbLine[46:54].strip()) - z_translate
                    new_x = x*rotationMatrix[0][0] + y*rotationMatrix[1][0] + z*rotationMatrix[2][0]
                    new_y = x*rotationMatrix[0][1] + y*rotationMatrix[1][1] + z*rotationMatrix[2][1]
                    new_z = x*rotationMatrix[0][2] + y*rotationMatrix[1][2] + z*rotationMatrix[2][2]
                    pdbLine = '%s%8.3f%8.3f%8.3f%s' %(pdbLine[0:30], new_x, new_y, new_z, pdbLine[54:len(pdbLine)])
                transformedPDBFile.write(pdbLine)
            pdbFile.close()
            transformedPDBFile.close()

    def rotationDictExtractor(self,phi, theta, psi):
        rotationMatrix = {}
        rotationMatrix[0] = [math.cos(theta)*math.cos(psi),math.cos(theta)*math.sin(psi),-math.sin(theta)]
        rotationMatrix[1] = [-math.cos(phi)*math.sin(psi) + math.sin(phi)*math.sin(theta)*math.cos(psi),math.cos(phi)*math.cos(psi) + math.sin(phi)*math.sin(theta)*math.sin(psi),math.sin(phi)*math.cos(theta)]
        rotationMatrix[2] = [math.sin(phi)*math.sin(psi) + math.cos(phi)*math.sin(theta)*math.cos(psi),-math.sin(phi)*math.cos(psi) + math.cos(phi)*math.sin(theta)*math.sin(psi),math.cos(phi)*math.cos(theta)]
        return rotationMatrix
        
    def transposeRotationDictExtractor(self,phi, theta, psi):
        rotationMatrix = {}
        rotationMatrix[0] = [math.cos(theta)*math.cos(psi),-math.cos(phi)*math.sin(psi) + math.sin(phi)*math.sin(theta)*math.cos(psi),math.sin(phi)*math.sin(psi) + math.cos(phi)*math.sin(theta)*math.cos(psi)]
        rotationMatrix[1] = [math.cos(theta)*math.sin(psi),math.cos(phi)*math.cos(psi) + math.sin(phi)*math.sin(theta)*math.sin(psi),-math.sin(phi)*math.cos(psi) + math.cos(phi)*math.sin(theta)*math.sin(psi)]
        rotationMatrix[2] = [-math.sin(theta),math.sin(phi)*math.cos(theta),math.cos(phi)*math.cos(theta)]
        return rotationMatrix
#########################################################################################################################
    def interfaceMatchCheck(self,dContact,matchDictLeft,matchDictRight):
        contact = 0
        for left in matchDictLeft.keys():
            for right in matchDictRight.keys():
                if dContact.has_key(left+right):
                    contact += 1
        if self.contact_Count <= contact:
            return True
        else:
            return False
    
    def matchThresholdCheck(self,partner):
        partnerName = partner[0]
        alignDict = partner[1]
        temp = partnerName.split("_")
        key = temp[0]+"_"+temp[1]
        #readhotspot file
        hotspotList = []
        hotsPath = self.hotspotPath + "/hotspot%s" % partnerName[:6]
        if os.path.exists(hotsPath):
            filehnd = open(hotsPath,"r")
            for line in filehnd.readlines():
                if not (line[0] == "#"):
                    hotspotList.append(line.strip().split()[0])
            filehnd.close()
        solutionList = []
        proteinSize = float(self.templateSize[key])
        if proteinSize <= 0:
            return -1
        for element in alignDict.keys():
            matchCount = alignDict[element][0]
            matchScore = (matchCount/proteinSize)*100
            hotspotAnalyze = self.hotspotAnalysis(alignDict[element][3],hotspotList) # returns 1 if interface successfully passed the hotspot test and 0 if it fails
            tmscore = alignDict[element][4]
            #print "Hotspot: %s , Match Count: %s , TM-score: %s\n" % (str(hotspotAnalyze), str(matchCount), str(tmscore))
            if hotspotAnalyze == 1 and matchCount >= self.MINIMUM_RESIDUE_MATCH_COUNT and tmscore >= self.TM_SCORE_THRESHOLD:
                if  proteinSize > self.template_Residue_Count:
                    if matchScore > (self.MINIMUM_RESIDUE_MATCH_PERCENTAGE-self.DIFF_PERCENTAGE):
                        solutionList.append(element)
                elif proteinSize <= self.template_Residue_Count:
                    if matchScore > self.MINIMUM_RESIDUE_MATCH_PERCENTAGE:
                        solutionList.append(element)
                else:
                    continue
        return solutionList,alignDict
    
    def alignmentDict(self,target,interface,chain):
        alignDict = -1

        #if target == "pdb1" or target == "pdb2": # when local alignment is donr for uploaded homolgy models
        fileName = "alignment/%s_%s_%s" % (interface,chain,target)
        if not (os.path.exists(fileName)):
            return alignDict
        fhn = open(fileName)
        alignDict = pickle.loads(fhn.read())
        fhn.close()
        # else:
        #     try:
        #         self.cur.execute("SELECT content FROM alignment_result where interface=%s && chain=%s && target=%s && method='TM-align'",(interface,chain,target))
        #         row = self.cur.fetchone()
        #         alignDict = pickle.loads(row[0])
        #     except:
        #         self.con.rollback() 
        return alignDict

    def hotspotAnalysis(self,matchDict,hotspotList):
        option = self.hotspotCriterion
        count = self.hotspotCount
        hotspotNum = 0
        if option == 0: #no hotspot needed
            return 1
        elif option == 1: #there should be count many hotspots in the matchDict
            for hotspot in hotspotList:
                if matchDict.has_key(hotspot):
                    hotspotNum += 1
            if count <= hotspotNum:
                return 1
            else:
                return 0
        elif option == 2: #there should be count many hotspots and corresponding residues should also be the same type
            for hotspot in hotspotList:
                if matchDict.has_key(hotspot) and matchDict[hotspot][2] == hotspot[2]: #same residue ex.A.S.164, we now took S
                    hotspotNum += 1
            if count <= hotspotNum:
                return 1
            else:
                return 0
        elif option == 3: # there should be count many hotspots and corresponding residues should be from the same class
            #classes
            #Hydrophobic - A,V,I,L,M,C = 0, Hydrophilic +charged, -charged, polar - K,R,H,D,E,S,T,P,N,Q = 1, Aromatic - F,Y,W = 2, Glycine - G = 3
            classes = {"A":0,"V":0,"I":0,"L":0,"M":0,"C":0,"K":1,"R":1,"H":1,"D":1,"E":1,"S":1,"T":1,"P":1,"N":1,"Q":1,"F":2,"Y":2,"W":2,"G":3}
            for hotspot in hotspotList:
                if matchDict.has_key(hotspot) and classes[matchDict[hotspot][2]] == classes[hotspot[2]]:
                    hotspotNum += 1
            if count <= hotspotNum:
                return 1
            else:
                return 0
        else:
            return 0 #actually hotspot criterian does not exist

