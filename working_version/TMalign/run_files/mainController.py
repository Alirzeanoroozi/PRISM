#calls all necessary classes automatically
#Written by Alper Baspinar
from pdbDownload import PDBdownload #pdbdownloader
from preProcessor import PreProcessor #preProcessor
from surfaceExtractor import SurfaceExtractor #surfaceExtractor
from structuralAlignmentTM import StructuralAligner #structuralAligner
from transformationFiltering import TransformFilter #transformFiltering
from flexibleRefinementRosetta import FlexibleRefinement #fiberdock step
from checkTemplate import TemplateChecker #checks if template exists and creates if needed
# from mysqlWriter import MysqlWriter #write to the database
# from databaseChecker import DatabaseChecker #check if the value exist or not saves time
import os

class Controller:
    def __init__(self,jobId):
        currentPath = os.getcwd()
        os.chdir("run_files")
        workPath = "../jobs/%s" % jobId
        filePath = "%s/results.php" % (workPath)
        listPath = "%s/lists" % (workPath)
        #prints indication of the module starts
        print("PDB download stage started...")
        leftTarget,rightTarget,templateList = PDBdownload(workPath).PDBdownloader()
        #checks template
        tempCheck = TemplateChecker(workPath,templateList).checker()
        checker = tempCheck[0]
        templateList = tempCheck[1]
        if checker == 1 or checker == 2:
            leftTarget = leftTarget[0:100] #at most 100 element will be evaluated for webserver
            rightTarget = rightTarget[0:100]
            prePdbList = leftTarget+rightTarget
            pdbList = []
            for p in prePdbList:
                if len(p) >= 4:
                    pdbList.append(p[0:4])
            pdbList = list(set(pdbList))
            print("PDB download stage finished...")
            print("PreProcess stage started...")
            try:
                print('Preparing proteins...')
                print('leftTarget,rightTarget,workPath:', leftTarget, rightTarget, workPath)
                leftTarget,rightTarget = PreProcessor(leftTarget,rightTarget,workPath).prepareProtein()
            except Exception as e:
                print e
                pdbList = []
            #a new module designed for network runs
            tableEntry = []
            if checker == 1:
                # leftTarget,rightTarget,tableEntry = DatabaseChecker(leftTarget,rightTarget).checker()
                pass
            pdbList = leftTarget+rightTarget
            pdbList = list(set(pdbList))
            print("PreProcess stage finished...")
            if checker == 1:
                # MysqlWriter(0,pdbList)
                pass
            print("SurfaceExtraction stage started...")
            try:
                pdbList = SurfaceExtractor(pdbList,workPath).surfaceExtractor()
            except Exception as e:
                print e
                pdbList = pdbList
            print("SurfaceExtraction stage finished...")
            #params.txt should be in workpath for multiprot to run correctly.
            if not (os.path.exists(workPath+"/params.txt")):
                os.system("cp %s %s/" % ("params.txt",workPath))
            print("Structural Alignment stage started...")
            StructuralAligner(pdbList,templateList,workPath)
            print("Structural Alignment stage finished...")
            print("Transformation Filtering stage started...")
            try:
                passedInterfaces = TransformFilter(leftTarget,rightTarget,templateList,workPath).transformer()
            except Exception as e:
                print e
                passedInterfaces = []
            print("Transformation Filtering stage finished...")
            #if checker == 1:
            #    MysqlWriter(2,passedInterfaces)    
            print("Flexible Refinement stage started...")
            energy_Structure = FlexibleRefinement(workPath,jobId).refiner()
            print("Flexible Refinement stage finished...")
            if checker == 1:
                # MysqlWriter(3,energy_Structure)
                # MysqlWriter(1,[leftTarget,rightTarget])
                pass
        elif checker == 0:
            print("Template Generation Failed...")
        # MysqlWriter(5,jobId)
        pass
        #remove folders in the work directory 
        #os.system("rm -r %s/*/" % (workPath))
        os.chdir(currentPath)
