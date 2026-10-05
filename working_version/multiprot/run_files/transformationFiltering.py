#!/usr/bin/env python
#Written by Alper Baspinar
import string,os,math,ConfigParser
import MySQLdb as mdb
import pickle
class TransformFilter:
    global leftTarget,rightTarget,templateList,MINIMUM_RESIDUE_MATCH_COUNT,MINIMUM_RESIDUE_MATCH_PERCENTAGE,MINIMUM_HOTSPOT_MATCH_NUMBER,DIFF_PERCENTAGE,multiprotOutPath,contactPath,hotspotPath,templateSize,multiprotCount,hotspotCriterion,hotspotCount,template_Residue_Count,contact_Count,clashing_distance,max_clashing_count,interfacePath,currentPath
    def __init__(self,leftTarget,rightTarget,templateList,workPath):
    	if os.path.exists("%s/fiberdock" % workPath):
    		return
	db_f = open("../config.inc", "r")
	db_f.readline()
	my_host = db_f.readline().split("'")[1]
	my_user = db_f.readline().split("'")[1]
	my_pass = db_f.readline().split("'")[1]
	my_db = db_f.readline().split("'")[1]
	db_f.close()
	self.con = mdb.connect(host=my_host, user=my_user, passwd=my_pass, db=my_db)
	self.cur = self.con.cursor()

	self.leftTarget = leftTarget
	self.rightTarget = rightTarget
        self.templateSize = {}
        self.templateList = templateList
        self.currentPath = os.getcwd()
        os.chdir(workPath)
        config = ConfigParser.ConfigParser() #reads configuration datas
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
        self.interfacePath = config.get('Structural_Alignment','interface_path')
        self.template_Residue_Count = config.getfloat('Transformation_Filtering','template_residue_count')
        self.contact_Count = config.getfloat('Transformation_Filtering','contact_count')
        self.clashing_distance = config.getfloat('Transformation_Filtering','clashing_distance')
        self.max_clashing_count = config.getfloat('Transformation_Filtering','max_clashing_count')
        if not (os.path.exists("transformation")):
		os.mkdir("transformation",0777)

    def transformer(self):
        filehnd = open("transformation/passedFiles","w") # if interface is successfull candidate write on the file
	#filehnd2 = open("transformation/passedInterfaces","w")
	passedInterfaces = []
        for interface in self.templateList:
            interface = interface[:6]
	    temp = self.filtering(interface)
            for passed in temp[0]:
                line = passed[0]+"\t"+passed[1]+"\n"
                filehnd.writelines(line)
	    for passInt in temp[1]:
		passedInterfaces.append(passInt)
		line = passInt+"\n"
		#filehnd2.writelines(line)
        filehnd.close()
	#filehnd2.close()
	os.chdir(self.currentPath)
	if self.con:
		self.con.close()
	return passedInterfaces

    def contactDict(self,interface):
        try:
            contactFile = open(self.contactPath + "/"+interface+'.txt','r')
        except:
            print "Interface %s does not exist..." % interface
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

            leftPartner = []
            rightPartner = []
            for index in range(len(self.leftTarget)):
            	left1 = "%s_%s_%s" % (interface,interface[4],self.leftTarget[index])
		multiDict1 = self.multiprotDict(self.leftTarget[index],interface,interface[4])
		left2 = "%s_%s_%s" % (interface,interface[4],self.rightTarget[index])
		multiDict2 = self.multiprotDict(self.rightTarget[index],interface,interface[4])
		right1 = "%s_%s_%s" % (interface,interface[5],self.rightTarget[index])
		multiDict3 = self.multiprotDict(self.rightTarget[index],interface,interface[5])
		right2 = "%s_%s_%s" % (interface,interface[5],self.leftTarget[index])
		multiDict4 = self.multiprotDict(self.leftTarget[index],interface,interface[5]) 
		
		if multiDict1 != -1 and multiDict3 != -1:
			check = True
			if len(self.leftTarget[index]) == 5 and len(self.rightTarget[index]) == 5:
				check = interface != (self.leftTarget[index]+self.rightTarget[index][-1])
			else:
				check = True
			if check:
				leftPartner.append([left1,multiDict1])
				rightPartner.append([right1,multiDict3])
		if multiDict2 != -1 and multiDict4 != -1:
			check = True
                        if len(self.leftTarget[index]) == 5 and len(self.rightTarget[index]) == 5:
				check = interface != (self.rightTarget[index]+self.leftTarget[index][-1])
                        else:
                                check = True
			if check:
				leftPartner.append([left2,multiDict2])
				rightPartner.append([right2,multiDict4])
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
		    solutionListLeft = leftMatchDict[index][0][0]
		    multiDictLeft = leftMatchDict[index][0][1]
		    left = leftMatchDict[index][1]
		    solutionListRight = rightMatchDict[index][0][0]
		    multiDictRight = rightMatchDict[index][0][1]
		    right = rightMatchDict[index][1]
		    for solLeft in solutionListLeft:
		  	for solRight in solutionListRight:
		            if self.interfaceMatchCheck(dContact,multiDictLeft[solLeft][3],multiDictRight[solRight][3]):
		                leftTarget = self.rotateTarget(left,solLeft,multiDictLeft[solLeft])
		                rightTarget = self.rotateTarget(right,solRight,multiDictRight[solRight])
		                if self.overlap(leftTarget,rightTarget):
		                    passedList.append([leftTarget,rightTarget])

            return passedList,passedInterfaces
    def overlap(self,leftTarget,rightTarget):
        clash_count = 0
        leftTarget = "transformation/"+leftTarget
        rightTarget = "transformation/"+rightTarget
        if not(os.path.exists(leftTarget)):
            print "Left Target %s does not exist." % leftTarget
            return False
        elif not(os.path.exists(rightTarget)):
            print "Right Target %s does not exist." % rightTarget
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
        phi, theta, psi, x_translate, y_translate, z_translate = transV
        if refMol == 0:
		    rotationMatrix = self.rotationDictExtractor(phi + math.pi, theta + math.pi, psi + math.pi)
		    pdbFile = open("preprocess/"+target, 'r')
		    transformedPDBFile = open("transformation/"+transformPath, 'w')
		    for pdbLine in pdbFile:
			    if pdbLine[0:4] == 'ATOM':
				    x = float(pdbLine[30:38].strip())
				    y = float(pdbLine[38:46].strip())
				    z = float(pdbLine[46:54].strip())
				    new_x = x*rotationMatrix[0][0] + y*rotationMatrix[1][0] + z*rotationMatrix[2][0] + x_translate
				    new_y = x*rotationMatrix[0][1] + y*rotationMatrix[1][1] + z*rotationMatrix[2][1] + y_translate
				    new_z = x*rotationMatrix[0][2] + y*rotationMatrix[1][2] + z*rotationMatrix[2][2] + z_translate
				    pdbLine = '%s%8.3f%8.3f%8.3f%s' %(pdbLine[0:30], new_x, new_y, new_z, pdbLine[54:len(pdbLine)])
                	    transformedPDBFile.write(pdbLine)
		    pdbFile.close()
		    transformedPDBFile.close()
        elif refMol == 1:
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
    	multiDict = partner[1]
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
        for element in multiDict.keys():
            matchCount = multiDict[element][0]
            matchScore = (matchCount/proteinSize)*100
            hotspotAnalyze = self.hotspotAnalysis(multiDict[element][3],hotspotList) # returns 1 if interface successfully passed the hotspot test and 0 if it fails
            if hotspotAnalyze == 1 and matchCount >= self.MINIMUM_RESIDUE_MATCH_COUNT:
                if  proteinSize > self.template_Residue_Count:
                    if matchScore > (self.MINIMUM_RESIDUE_MATCH_PERCENTAGE-self.DIFF_PERCENTAGE):
                        solutionList.append(element)
                elif proteinSize <= self.template_Residue_Count:
                    if matchScore > self.MINIMUM_RESIDUE_MATCH_PERCENTAGE:
                        solutionList.append(element)
                else:
                    continue
        return solutionList,multiDict
    
    def multiprotDict(self,target,interface,chain):
       multiDict = -1
       if target == "pdb1" or target == "pdb2":
		fileName = "alignment/%s_%s_%s" % (interface,chain,target)
		if not (os.path.exists(fileName)):
			return multiDict
		fhn = open(fileName)
		multiDict = pickle.loads(fhn.read())
		fhn.close()
       else:
		try:
			self.cur.execute("SELECT content FROM multiprot where interface=%s && chain=%s && target=%s",(interface,chain,target))
			row = self.cur.fetchone()
			multiDict = pickle.loads(row[0])
		except:
			self.con.rollback() 
       return multiDict

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

