#calls all necessary classes automatically
#Written by Alper Baspinar
from pdbDownload import PDBdownload #pdbdownloader
from preProcessor import PreProcessor #preProcessor
from surfaceExtractor import SurfaceExtractor #surfaceExtractor
from structuralAlignment import StructuralAligner #structuralAligner
from transformationFiltering import TransformFilter #transformFiltering
from flexibleRefinement import FlexibleRefinement #fiberdock step
from checkTemplate import TemplateChecker #checks if template exists and creates if needed
from htmlWriter import HtmlWriter #create htmlfile
from mysqlWriter import MysqlWriter #write to the database
from databaseChecker import DatabaseChecker #check if the value exist or not saves time
from sendMail import MailSender # check if user provided a mail address and sends mail if there is one
import os,sys

class Controller:
	def __init__(self,jobId):
    		currentPath = os.getcwd()
    		os.chdir("run_files")
    		
    		workPath = "../jobs/%s" % jobId
		filePath = "%s/results.php" % (workPath)
		listPath = "%s/lists"	% (workPath)
		#prints indication of the module starts
		HtmlWriter(filePath,listPath,-1,[])
		print("PDB download stage started...")
		HtmlWriter(filePath,listPath,0,[])
		leftTarget,rightTarget,templateList = PDBdownload(workPath).PDBdownloader()
		#checks template
		tempCheck = TemplateChecker(workPath,templateList).checker()
		checker = tempCheck[0]
		templateList = tempCheck[1]
		if checker == 1 or checker == 2:
			leftTarget = leftTarget[0:10] #at most 10 element will be evaluated for webserver
			rightTarget = rightTarget[0:10]
			prePdbList = leftTarget+rightTarget
			pdbList = []
			for p in prePdbList:
				if len(p) >= 4:
					pdbList.append(p[0:4])	
			pdbList = list(set(pdbList))
			print("PDB download stage finished...")
			HtmlWriter(filePath,listPath,1,pdbList)
			print("PreProcess stage started...")
			HtmlWriter(filePath,listPath,2,prePdbList)
			try:
				leftTarget,rightTarget = PreProcessor(leftTarget,rightTarget,workPath).prepareProtein()
			except Exception as e:
				print e
				pdbList = []
			#a new module designed for network runs
			tableEntry = []
			if checker == 1:
				leftTarget,rightTarget,tableEntry = DatabaseChecker(leftTarget,rightTarget).checker()
			pdbList = leftTarget+rightTarget
			pdbList = list(set(pdbList))
			print("PreProcess stage finished...")
			HtmlWriter(filePath,listPath,3,pdbList)
			if checker == 1:
				MysqlWriter(0,pdbList)
			print("SurfaceExtraction stage started...")
			HtmlWriter(filePath,listPath,4,pdbList)
			try:
				pdbList = SurfaceExtractor(pdbList,workPath).surfaceExtractor()
			except Exception as e:
				print e
				pdbList = pdbList
			print("SurfaceExtraction stage finished...")
			HtmlWriter(filePath,listPath,5,pdbList)
			#params.txt should be in workpath for multiprot to run correctly.
			if not (os.path.exists(workPath+"/params.txt")):
				os.system("cp %s %s/" % ("params.txt",workPath))
		
			print("Structural Alignment stage started...")
			HtmlWriter(filePath,listPath,6,pdbList)
			StructuralAligner(pdbList,templateList,workPath)
			print("Structural Alignment stage finished...")
			HtmlWriter(filePath,listPath,7,pdbList)
		
			print("Transformation Filtering stage started...")
			HtmlWriter(filePath,listPath,8,pdbList)
			try:
				passedInterfaces = TransformFilter(leftTarget,rightTarget,templateList,workPath).transformer()
			except Exception as e:
				print e
				passedInterfaces = []
			print("Transformation Filtering stage finished...")
			HtmlWriter(filePath,listPath,9,pdbList)
			#if checker == 1:
			#	MysqlWriter(2,passedInterfaces)	
			print("Flexible Refinement stage started...")
			HtmlWriter(filePath,listPath,10,pdbList)
			energy_Structure = FlexibleRefinement(workPath,jobId).refiner()
			print("Flexible Refinement stage finished...")
			HtmlWriter(filePath,listPath,11,pdbList)
			HtmlWriter(filePath,listPath,12,[tableEntry,energy_Structure])
			if checker == 1:
				MysqlWriter(3,energy_Structure)
				MysqlWriter(1,[leftTarget,rightTarget])
		elif checker == 0:
			HtmlWriter(filePath,listPath,12,[[],[]])
			print("Template Generation Failed...")
		MysqlWriter(5,jobId)
		MailSender(jobId)
		#remove folders in the work directory 
		os.system("rm -r %s/*/" % (workPath))
		os.chdir(currentPath)
