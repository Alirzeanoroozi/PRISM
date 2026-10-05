#!/usr/bin/env python
#Written by Alper Baspinar
import os,ConfigParser

from ftplib import FTP

#downloads protein pdbs from the pdb databank
class PDBdownload:
    #constructor: workPath should be provided.
    global currentPath,pdbPath,workPath
    def __init__(self,workPath):
    	self.currentPath = os.getcwd() #holds the current directory.
    	self.workPath = os.path.abspath(workPath)
    	os.chdir(self.workPath) #changes directory to workPath       
    	config = ConfigParser.ConfigParser()
    	config.read('prism.ini') #use configParser to get configuration info. from 'prism.ini'
    	self.pdbPath = os.path.abspath(config.get('Pdb_Folder','pdb_path')) #defines where to download pdbs wrt workpath and changes into abs
    	     
    def PDBdownloader(self):
	preLeftTarget = []
	preRightTarget = []
    	preTemplateList = []
    	#checks whether pairList and templateList exists
    	pairListPath = "lists/pair_list"
    	templateListPath = "lists/template_list"
    	if  not (os.path.exists(pairListPath)) or not (os.path.exists(templateListPath)):
    		print "pdblist or templatelist does not exists..."
    		os.chdir(self.currentPath)
		return preLeftTarget,preRightTarget,preTemplateList
    	else:
    		filehnd = open(pairListPath,"r")
    		for line in filehnd.readlines():
			line = line.strip()
			line = line.split()
			if(len(line) == 2):
				if len(line[0]) >= 4 and len(line[1]) >= 4:
					preLeftTarget.append(line[0])
					preRightTarget.append(line[1])
    		filehnd.close()
    		filehnd = open(templateListPath,"r")
    		for line in filehnd.readlines():
    			preTemplateList.append(line.strip()[:6])
    		filehnd.close()
	#prechecks preLeftTarget and preRightTarget before process
	leftTarget = []
	for pdb in preLeftTarget:
                newPdb = pdb[:4].lower() #the first 4 letter of a pdb should be lower case 
                extra = sorted(set(pdb[4:])) #pdb chains must be unique and in ascending order
                for e in extra:
			newPdb += e
                leftTarget.append(newPdb)

        rightTarget = []
	for pdb in preRightTarget:
                newPdb = pdb[:4].lower() #the first 4 letter of a pdb should be lower case 
                extra = sorted(set(pdb[4:])) #pdb chains must be unique and in ascending order
                for e in extra:
                    	newPdb += e
                rightTarget.append(newPdb)

        if os.path.exists("preprocess"):
    		os.chdir(self.currentPath)
    		return leftTarget,rightTarget,preTemplateList 
    			
        if len(leftTarget+rightTarget) == 0:
        	print "Pdblist does not contain any appropriate pdb"
        	os.chdir(self.currentPath)
        	return leftTarget,rightTarget,preTemplateList
       	
	#remove duplicate entries
	l = []
	r = []
	tt = []
	for i in range(len(leftTarget)):
		t = str(leftTarget[i])+str(rightTarget[i])
		tr = str(rightTarget[i])+str(leftTarget[i])
		if not (t in tt):
			tt.append(t)
			tt.append(tr)
			l.append(leftTarget[i])
			r.append(rightTarget[i])
			
	leftTarget = l
	rightTarget = r
	
        localPdb = "%s/pdb" % self.workPath
        if not (os.path.exists(localPdb)):
        	os.mkdir(localPdb,0777)	
        #first checks if the pdb already exist and then downloads
	leftList = []
	rightList = []
	for index in range(len(leftTarget)):
		proteinNameLeft = leftTarget[index]
		proteinNameRight = rightTarget[index]
		check1 = False
		check2 = False
		if os.path.exists("%s/%s.pdb" % (localPdb,proteinNameLeft[0:4])):
			check1 = True
		elif os.path.exists("%s/%s.pdb" % (self.pdbPath,proteinNameLeft[0:4])):#pdbList might contain protein with chains, but pdb databank needs length 4 inputs
			check1 = True
			os.system("cp %s/%s.pdb %s/%s.pdb" % (self.pdbPath,proteinNameLeft[0:4],localPdb,proteinNameLeft[0:4]))
		else:
			os.chdir(self.pdbPath) # changes directory to where pdb files are located wrt workPath.
			#try and except to fetch pdb files each download creates a new ftp connection for fail tolerance
		        ftp = FTP('ftp.wwpdb.org')
		        try:
		            #initialize ftp connection to download pdb files
		            ftp.login()
		            ftp.cwd('pub/pdb/data/structures/all/pdb')
		            ftp.retrbinary('RETR pdb'+proteinNameLeft[0:4]+'.ent.gz', open('pdb'+proteinNameLeft[0:4]+'.ent.gz', 'wb').write)
		            os.system("gunzip pdb%s.ent.gz" % proteinNameLeft[0:4])#gunzip pdbs because it is in .ent.gz format
		            os.system("mv pdb%s.ent %s.pdb" % (proteinNameLeft[0:4],proteinNameLeft[0:4]))#change name from pdb____.ent to ____.pdb
		            check1 = True
		            os.system("cp %s/%s.pdb %s/%s.pdb" % (self.pdbPath,proteinNameLeft[0:4],localPdb,proteinNameLeft[0:4]))
		        except:
		            print "could not fetch the %s from online pdb database" % proteinNameLeft[0:4]
		            os.system("rm *.ent.gz")
		        ftp.quit()
		
		if os.path.exists("%s/%s.pdb" % (localPdb,proteinNameRight[0:4])):
			check2 = True
		elif os.path.exists("%s/%s.pdb" % (self.pdbPath,proteinNameRight[0:4])):#pdbList might contain protein with chains, but pdb databank needs length 4 inputs
			check2 = True
			os.system("cp %s/%s.pdb %s/%s.pdb" % (self.pdbPath,proteinNameRight[0:4],localPdb,proteinNameRight[0:4]))
		else:
			if os.getcwd() != self.pdbPath: 
				os.chdir(self.pdbPath) # changes directory to where pdb files are located wrt workPath.
			#try and except to fetch pdb files each download creates a new ftp connection for fail tolerance
		        ftp = FTP('ftp.wwpdb.org')
		        try:
		            #initialize ftp connection to download pdb files
		            ftp.login()
		            ftp.cwd('pub/pdb/data/structures/all/pdb')
		            ftp.retrbinary('RETR pdb'+proteinNameRight[0:4]+'.ent.gz', open('pdb'+proteinNameRight[0:4]+'.ent.gz', 'wb').write)
		            os.system("gunzip pdb%s.ent.gz" % proteinNameRight[0:4])#gunzip pdbs because it is in .ent.gz format
		            os.system("mv pdb%s.ent %s.pdb" % (proteinNameRight[0:4],proteinNameRight[0:4]))#change name from pdb____.ent to ____.pdb
		            check2 = True
		            os.system("cp %s/%s.pdb %s/%s.pdb" % (self.pdbPath,proteinNameRight[0:4],localPdb,proteinNameRight[0:4]))
		        except:
		            print "could not fetch the %s from online pdb database" % proteinNameRight[0:4]
		            os.system("rm *.ent.gz")
		        ftp.quit()
                if check1 and check2:
                	leftList.append(proteinNameLeft)
                	rightList.append(proteinNameRight)
                	
        os.chdir(self.currentPath)
        return leftList,rightList,preTemplateList

