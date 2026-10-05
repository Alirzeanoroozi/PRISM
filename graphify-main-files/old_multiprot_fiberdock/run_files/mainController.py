#calls all necessary classes automatically
#Written by Alper Baspinar
#CLI/fiberdock port: keeps multiprot alignment + fiberdock refinement;
#MySQL/mail/html-progress calls disabled for command-line runs.
from pdbDownload import PDBdownload #pdbdownloader
from preProcessor import PreProcessor #preProcessor
from surfaceExtractor import SurfaceExtractor #surfaceExtractor
from structuralAlignment import StructuralAligner #structuralAligner (multiprot)
from transformationFiltering import TransformFilter #transformFiltering
from flexibleRefinement import FlexibleRefinement #fiberdock refinement step
from checkTemplate import TemplateChecker #checks if template exists and creates if needed
# from htmlWriter import HtmlWriter #webserver progress page (disabled for CLI)
# from mysqlWriter import MysqlWriter #write to the database (disabled for CLI)
# from databaseChecker import DatabaseChecker #db dedup (disabled for CLI)
# from sendMail import MailSender #completion mail (disabled for CLI)
import os

class Controller:
    def __init__(self,jobId):
        currentPath = os.getcwd()
        os.chdir("run_files")

        workPath = "../jobs/%s" % jobId
        listPath = "%s/lists" % (workPath)

        print("PDB download stage started...")
        leftTarget,rightTarget,templateList = PDBdownload(workPath).PDBdownloader()
        #checks template (generates missing templates via templateGenerator)
        tempCheck = TemplateChecker(workPath,templateList).checker()
        checker = tempCheck[0]
        templateList = tempCheck[1]
        if checker == 1 or checker == 2:
            leftTarget = leftTarget[0:100] #webserver capped at 10; lifted to 100 for CLI
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
                leftTarget,rightTarget = PreProcessor(leftTarget,rightTarget,workPath).prepareProtein()
            except Exception as e:
                print e
                pdbList = []
            #database dedup step (network runs only) -- disabled for CLI
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

            print("Flexible Refinement stage started...")
            energy_Structure = FlexibleRefinement(workPath,jobId).refiner()
            print("Flexible Refinement stage finished...")
            print("Results (energy / structure):")
            for row in energy_Structure:
                print(row)
            if checker == 1:
                # MysqlWriter(3,energy_Structure)
                # MysqlWriter(1,[leftTarget,rightTarget])
                pass
        elif checker == 0:
            print("Template Generation Failed...")
        # MysqlWriter(5,jobId)
        # MailSender(jobId)
        #keep intermediate folders for CLI inspection (webserver cleaned them up)
        #os.system("rm -r %s/*/" % (workPath))
        os.chdir(currentPath)
