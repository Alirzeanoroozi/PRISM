import os,math
import configparser
import pickle
from Bio.PDB import PDBParser


os.makedirs("transformation", exist_ok=True)

class TransformFilter:
    def __init__(self):
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
        self.template_Residue_Count = config.getfloat('Transformation_Filtering','template_residue_count')
        self.contact_Count = config.getfloat('Transformation_Filtering','contact_count')
        self.clashing_distance = config.getfloat('Transformation_Filtering','clashing_distance')
        self.max_clashing_count = config.getfloat('Transformation_Filtering','max_clashing_count')
        self.TM_SCORE_THRESHOLD = config.getfloat('Structural_Alignment','tm_score_threshold')

    def transformer(self, queries, templates):
        for interface in templates:
            temp = self.filtering(interface)

    def filtering(self, interface):
        sizePathLeft =  f"template/interfaces/{interface}_{interface[4]}.int"
        sizePathRight = f"template/interfaces/{interface}_{interface[5]}.int"

        contact_dict = self.contact_dict(interface)

        with open(sizePathLeft, "r") as f:
            left_key = f"{interface}_{interface[4]}"
            self.templateSize[left_key] = len(f.readlines())
        with open(sizePathRight, "r") as f:
            right_key = f"{interface}_{interface[5]}"
            self.templateSize[right_key] = len(f.readlines())

        leftPartner = []
        rightPartner = []
        for index in range(len(self.leftTarget)):
            left1 = f"{interface}_{interface[4]}_{self.leftTarget[index]}"
            align_dict1 = alignment_dict(self.leftTarget[index], interface, interface[4])
            left2 = f"{interface}_{interface[4]}_{self.rightTarget[index]}"
            align_dict2 = alignment_dict(self.rightTarget[index], interface, interface[4])
            right1 = f"{interface}_{interface[5]}_{self.rightTarget[index]}"
            align_dict3 = alignment_dict(self.rightTarget[index], interface, interface[5])
            right2 = f"{interface}_{interface[5]}_{self.leftTarget[index]}"
            align_dict4 = alignment_dict(self.leftTarget[index], interface, interface[5]) 

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
                    if self.interfaceMatchCheck(contact_dict,alignDictLeft[solLeft][3],alignDictRight[solRight][3]):
                        leftTarget = self.rotateTarget(left,solLeft,alignDictLeft[solLeft])
                        rightTarget = self.rotateTarget(right,solRight,alignDictRight[solRight])
                        if self.number_of_clashing_residues(leftTarget, rightTarget):
                            passedList.append([leftTarget,rightTarget])
        return passedList,passedInterfaces
        
    def rotateTarget(self, partner, sol,multiList):
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

#########################################################################################################################
    def interface_match_check(self, contact_dict, matchDictLeft, matchDictRight):
        contact = 0
        for left in matchDictLeft.keys():
            for right in matchDictRight.keys():
                if dContact.has_key(left+right):
                    contact += 1
        if self.contact_Count <= contact:
            return True
        else:
            return False
    
    def match_threshold_check(self, partner):
        partnerName = partner[0]
        align_results = partner[1]
        match_count, translation, rotationMat, match_dict, tmscore = align_results
        temp = partnerName.split("_")
        key = temp[0]+"_"+temp[1]
        #readhotspot file
        hotspotList = []
        hotsPath = f"template/hotspot/hotspot{partnerName[:6]}"
        if os.path.exists(hotsPath):
            with open(hotsPath, "r") as f:
                for line in f.readlines():
                    if not (line[0] == "#"):
                        hotspotList.append(line.strip().split()[0])
        solutionList = []
        proteinSize = float(self.templateSize[key])
        if proteinSize <= 0:
            return -1

        match_score = (match_count / proteinSize) * 100
        hotspotAnalyze = self.hotspot_analysis(match_dict, hotspotList) # returns 1 if interface successfully passed the hotspot test and 0 if it fails
        if hotspotAnalyze == 1 and match_count >= self.MINIMUM_RESIDUE_MATCH_COUNT and tmscore >= self.TM_SCORE_THRESHOLD:
            if  proteinSize > self.template_Residue_Count:
                if match_score > (self.MINIMUM_RESIDUE_MATCH_PERCENTAGE-self.DIFF_PERCENTAGE):
                    solutionList.append(partnerName)
            elif proteinSize <= self.template_Residue_Count:
                if match_score > self.MINIMUM_RESIDUE_MATCH_PERCENTAGE:
                    solutionList.append(partnerName)
        return solutionList, align_results
    
    def hotspot_analysis(self, match_dict, hotspot_list):
        option = self.hotspotCriterion
        count = self.hotspotCount
        hotspot_num = 0
        if option == 0: #no hotspot needed
            return 1
        elif option == 1: #there should be count many hotspots in the matchDict
            for hotspot in hotspot_list:
                if hotspot in match_dict:
                    hotspot_num += 1
            if count <= hotspot_num:
                return 1
            else:
                return 0
        elif option == 2: #there should be count many hotspots and corresponding residues should also be the same type
            for hotspot in hotspot_list:
                if hotspot in match_dict and match_dict[hotspot][2] == hotspot[2]: #same residue ex.A.S.164, we now took S
                    hotspot_num += 1
            if count <= hotspot_num:
                return 1
            else:
                return 0
        elif option == 3: # there should be count many hotspots and corresponding residues should be from the same class
            #classes
            #Hydrophobic - A,V,I,L,M,C = 0, Hydrophilic +charged, -charged, polar - K,R,H,D,E,S,T,P,N,Q = 1, Aromatic - F,Y,W = 2, Glycine - G = 3
            classes = {"A":0,"V":0,"I":0,"L":0,"M":0,"C":0,"K":1,"R":1,"H":1,"D":1,"E":1,"S":1,"T":1,"P":1,"N":1,"Q":1,"F":2,"Y":2,"W":2,"G":3}
            for hotspot in hotspot_list:
                if hotspot in match_dict and classes[match_dict[hotspot][2]] == classes[hotspot[2]]:
                    hotspot_num += 1
            if count <= hotspot_num:
                return 1
            else:
                return 0
        else:
            return 0 #actually hotspot criterian does not exist

def number_of_clashing_residues(leftTarget, rightTarget):
    clashing_distance = 3
    max_clashing_count = 5
    clash_count = 0

    left_ca_coordinates = read_target("transformation/"+leftTarget)
    right_ca_coordinates = read_target("transformation/"+rightTarget)
    for left_ca_coordinate in left_ca_coordinates:
        for right_ca_coordinate in right_ca_coordinates:
            if get_distance(left_ca_coordinate, right_ca_coordinate) < clashing_distance:
                clash_count += 1
                if max_clashing_count <= clash_count:
                    return False
    return True

def read_target(path):
    parser = PDBParser(QUIET=True)
    ca_coordinates = []
    structure = parser.get_structure("target", path)
    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    ca_atom = residue['CA']
                    coord = list(ca_atom.get_coord())
                    ca_coordinates.append(coord)
    return ca_coordinates

def get_distance(coord1, coord2):
    x1,y1,z1 = coord1
    x2,y2,z2 = coord2
    return ((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)**0.5

def contact_dict(interface):
    with open(f"template/contact/{interface}.txt", 'r') as f:
        return {line.strip().split()[0] + line.strip().split()[1]: 1 for line in f.readlines()}

def transpose_rotation_dict_extractor(phi, theta, psi):
    return {
        0: [math.cos(theta)*math.cos(psi),-math.cos(phi)*math.sin(psi) + math.sin(phi)*math.sin(theta)*math.cos(psi),math.sin(phi)*math.sin(psi) + math.cos(phi)*math.sin(theta)*math.cos(psi)],
        1: [math.cos(theta)*math.sin(psi),math.cos(phi)*math.cos(psi) + math.sin(phi)*math.sin(theta)*math.sin(psi),-math.sin(phi)*math.cos(psi) + math.cos(phi)*math.sin(theta)*math.sin(psi)],
        2: [-math.sin(theta),math.sin(phi)*math.cos(theta),math.cos(phi)*math.cos(theta)]
    }

def rotation_dict_extractor(phi, theta, psi):
    return {
        0: [math.cos(theta)*math.cos(psi),math.cos(theta)*math.sin(psi),-math.sin(theta)],
        1: [-math.cos(phi)*math.sin(psi) + math.sin(phi)*math.sin(theta)*math.cos(psi),math.cos(phi)*math.cos(psi) + math.sin(phi)*math.sin(theta)*math.sin(psi),math.sin(phi)*math.cos(theta)],
        2: [math.sin(phi)*math.sin(psi) + math.cos(phi)*math.sin(theta)*math.cos(psi),-math.sin(phi)*math.cos(psi) + math.cos(phi)*math.sin(theta)*math.sin(psi),math.cos(phi)*math.cos(theta)]
    }

def alignment_dict(target, interface, chain):
    with open(f"alignment/{interface}_{chain}_{target}.pkl", "r") as f:
        return pickle.load(f)
