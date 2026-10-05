#!/usr/bin/env python
#Written by Alper Baspinar
import os,sys,ConfigParser
import fiberdockInterfaceExtractor as fb
class FlexibleRefinement:
    global workPath,jobId,currentPath,hydrogenAdder,buildParams,fiberdock,nma,normal_modes,sizes,ca_Num,aDict
    def __init__(self,workPath,jobId):
        self.workPath = workPath
        self.jobId = jobId
        self.sizes = {}
        self.ca_Num = {}
        self.lastResidue = {}
        self.aDict = {}
        self.currentPath = os.getcwd()
        os.chdir(workPath)
        config = ConfigParser.ConfigParser() #reads configuration datas
        config.read('prism.ini')
        self.hydrogenAdder = config.get('External_Tools','addhydrogen')
        self.buildParams = config.get('External_Tools','fiberdockparams')
        self.fiberdock = config.get('External_Tools','fiberdock')
        self.nma = config.get('External_Tools','nma')
        self.normal_modes = config.getint('Flexible_Refinement','normal_modes')
        if not (os.path.exists("fiberdock")):
            os.mkdir("fiberdock",0777)
        if not (os.path.exists("fiberdock/energies")):
            os.mkdir("fiberdock/energies",0777) 
        if not (os.path.exists("fiberdock/structures")):
            os.mkdir("fiberdock/structures",0777)  

    def refiner(self):
        energy_Structure = []
        energyFile = "fiberdock/energies/fiberdock_energies"
        filehndOut = open(energyFile,"w")
        if not (os.path.exists("transformation/passedFiles")):
            print "passedFiles does not exists check other steps."
        else:
            hnd = open("fiberdock/zero-transformation","w")
            hnd.writelines("1 0 0 0 0 0 0")
            hnd.close()
            filehnd = open("transformation/passedFiles","r")
            for passed in filehnd.readlines():
                passed = passed.split()
                passed0 = passed[0].strip()
                passed1 = passed[1].strip()
                energy,structure = self.calculateEnergy(passed0,passed1)
                if energy != "-":
                    energy_Structure.append(["%s\t%s\t%s" % (passed0,passed1,energy),structure])
                    filehndOut.writelines("%s\t%s\t%s\n" % (passed0,passed1,energy))
            filehnd.close()
            filehndOut.close()
        if os.path.exists("fiberdock.log"):
            os.system("mv fiberdock.log fiberdock")
        os.chdir(self.currentPath)
        return energy_Structure

    def calculateEnergy(self,passed0,passed1):
        caPath1 = "fiberdock/%s.ca.pdb" % (passed0)
        caPath2 = "fiberdock/%s.ca.pdb" % (passed1)
        hbPath1 = "fiberdock/%s.HB" % (passed0)
        hbPath2 = "fiberdock/%s.HB" % (passed1)
        nmaPath1 = "fiberdock/%s.ca.nma" % (passed0)
        nmaPath2 = "fiberdock/%s.ca.nma" % (passed1)
 
        if not (os.path.exists(hbPath1)):
            hbOut = "fiberdock/%s.reduce" % passed0
            os.system("%s %s > %s" % (self.hydrogenAdder,"transformation/"+passed0,hbOut))
            if os.path.exists("transformation/%s.HB" % passed0):
                os.system("mv transformation/%s.HB %s" % (passed0,hbPath1))
        if not (os.path.exists(hbPath2)):
            hbOut = "fiberdock/%s.reduce" % passed1
            os.system("%s %s > %s" % (self.hydrogenAdder,"transformation/"+passed1,hbOut))
            if os.path.exists("transformation/%s.HB" % passed1):
                os.system("mv transformation/%s.HB %s" % (passed1,hbPath2))
        if not (os.path.exists(caPath1)) or not(self.sizes.has_key(passed0)):
            temp = self.createCaAtoms(passed0,caPath1)
            self.sizes[passed0] = temp[0]
            self.ca_Num[passed0] = temp[1]
        if not (os.path.exists(caPath2)) or not(self.sizes.has_key(passed1)):
            temp = self.createCaAtoms(passed1,caPath2)
            self.sizes[passed1] = temp[0]
            self.ca_Num[passed1] = temp[1]

        if not (os.path.exists(nmaPath1)):
            nmaOut = "fiberdock/%s.NMA" % passed0
            os.system("%s %s %s %d 3 10 > %s" % (self.nma,caPath1,nmaPath1,self.normal_modes,nmaOut))
            if os.path.exists("%s.ca.nma" % passed0):
                os.system("mv %s.ca.nma %s" % (passed0,nmaPath1))

        if not (os.path.exists(nmaPath2)):
            nmaOut = "fiberdock/%s.NMA" % passed1
            os.system("%s %s %s %d 3 10 > %s" % (self.nma,caPath2,nmaPath2,self.normal_modes,nmaOut))
            if os.path.exists("%s.ca.nma" % passed1):
            	os.system("mv %s.ca.nma %s" % (passed1,nmaPath2))
        
        energy = "-"
        structure = ""
        if self.sizes[passed0] >= self.sizes[passed1]:
            energy,structure = self.runFiberdock(passed0,passed1,0)
        else:
            energy,structure = self.runFiberdock(passed1,passed0,1)
        return str(energy),str(structure)

    def runFiberdock(self,receptor,ligand,ref):
        #build parameter zero-trial file first
        fiberdockOut = "fiberdock/%s_%s.fib" % (receptor,ligand) 
        os.system("%s fiberdock/%s.HB fiberdock/%s.HB U U Default fiberdock/zero-transformation fiberdock/%s_%s 0 50 0.80 1 glpk zero-trial 0.05 fiberdock/%s.ca.pdb fiberdock/%s.ca.nma fiberdock/%s.ca.pdb fiberdock/%s.ca.nma" % (self.buildParams,receptor,ligand,receptor,ligand,receptor,receptor,ligand,ligand))
        os.system("%s zero-trial > %s" % (self.fiberdock,fiberdockOut))
        #energy part
        fiberdockEnergy = "fiberdock/%s_%s.ref" % (receptor,ligand)
        energy = "-"
        if os.path.exists(fiberdockEnergy):
            filehnd = open(fiberdockEnergy,"r")
            for line in filehnd.readlines():
                if line.find("|") == -1:
                    continue
                else:
                    line = line.split("|")[1].strip()
                    if line == "glob":
                        continue
                    else:
                        try:
                            energy = float(line)
                            if energy > 0:
                                energy = "-"
                        except:
                            energy = "-"
            filehnd.close()
        fiberdockStructure = "fiberdock/%s_%s_1.ref.pdb" % (receptor,ligand)
        structureList = {0:[],1:[]}
        structureOut = "-"
        global_fib = "-"
        global_intR = "-"
        if os.path.exists(fiberdockStructure) and energy != "-":
            filehnd = open(fiberdockStructure,"r")
            i = 0
            ca_count = self.ca_Num[receptor]
            while i < ca_count:
                line = filehnd.readline()
                if line == "":
                    break
                structureList[0].append(line)
                atom = line[12:16].strip()
                if atom == "CA":
                    i += 1
            line = filehnd.readline()
            while  not (line == ""):
                atom = line[12:16].strip()
                if atom == "N":
                    break
                else:
                    structureList[0].append(line)
                line = filehnd.readline()
            while not (line == ""):
                structureList[1].append(line)
                line = filehnd.readline()
            filehnd.close()
            global_fib = "fiberdock_output/%s/" % (self.jobId)
            global_intR = "fiberdock_output/%s/" % (self.jobId)
            if not (os.path.exists("../../%s" % (global_fib))):
                os.mkdir("../../%s" % (global_fib),0777)
                os.system("chmod 775 ../../%s" % (global_fib))
            #parse name of receptor and ligand
            temp1 = receptor.split("_")
            temp2 = ligand.split("_")
            key1 = temp1[0]+","+temp1[2].split(".")[0]+","+temp2[2].split(".")[0]+","+str(energy)
            key2 = temp2[0]+","+temp2[2].split(".")[0]+","+temp1[2].split(".")[0]+","+str(energy)
            check = True
            if ref == 0:
                if self.aDict.has_key(key1):
                    check = False
                else:
                    self.aDict[key1] = 1
                outputName = temp1[0]+"_"+temp1[2].split(".")[0]+"_"+temp1[3]+"_"+temp2[2].split(".")[0]+"_"+temp2[3]
                structureOut = "fiberdock/structures/%s.fiberdock.pdb" % (outputName)
                global_fib += "%s.fiberdock.pdb" % (outputName) 
                global_intR += "%s.intRes.txt" % (outputName)
            else:
                if self.aDict.has_key(key2):
                    check = False
                else:
                    self.aDict[key2] = 1
                outputName = temp2[0]+"_"+temp2[2].split(".")[0]+"_"+temp2[3]+"_"+temp1[2].split(".")[0]+"_"+temp1[3]
                structureOut = "fiberdock/structures/%s.fiberdock.pdb" % (outputName)
                global_fib += "%s.fiberdock.pdb" % (outputName)
                global_intR += "%s.intRes.txt" % (outputName)
        
            if check:
                try:
                    fb.fiberdockInterfaceExtractor("../../"+global_fib,"../../"+global_intR,structureList[ref],structureList[1-ref])
                except:
                    print "Fiberdock output did not created."
            else:
                energy = "-"
    
        return energy,"%s" % (global_fib)

    def createCaAtoms(self,passed,caPath):
        if not(os.path.exists("transformation/"+passed)):
            return 0,0
        filehnd = open("transformation/"+passed,"r")
        cahnd = open(caPath,"w")
        size = 0
        caNum = 0
        for line in filehnd.readlines():
            if line[:3] == "END":
                break
            elif line[:4] == "ATOM":
                size += 1
                if line[12:16].strip() == "CA":
                    caNum += 1
                    cahnd.writelines(line)
        cahnd.writelines("END")
        cahnd.close()
        filehnd.close()
        return size,caNum


