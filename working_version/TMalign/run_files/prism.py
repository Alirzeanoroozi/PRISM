import os,sys,shutil
run_files_path = os.path.abspath(os.path.dirname(__file__))
project_root = os.path.dirname(run_files_path)
if not run_files_path in sys.path:
    sys.path.insert(1, run_files_path)
from mainController import Controller

if len(sys.argv) == 4:
    os.chdir(project_root)
    pdbListPath = sys.argv[1]
    templateListPath = sys.argv[2]
    jobId = sys.argv[3]
    currentPath = os.getcwd()
    #creates a folder for the workPath
    workPath = "jobs/%s" % jobId
    if os.path.exists(workPath) == False:
	os.makedirs(workPath,0777)
    #creates a folder for the lists
    listPath = workPath+"/lists"
    if os.path.exists(listPath) == False:
        os.mkdir(listPath,0777)
    #moves pdbList to workFolder
    pdbPath = listPath+"/pair_list"
    if not (os.path.exists(pdbPath)):
        shutil.copy2(pdbListPath,pdbPath)
    #moves templateList to workFolder
    templatePath = listPath+"/template_list"
    if not (os.path.exists(templatePath)):
        shutil.copy2(templateListPath,templatePath)
    
    configPath = "%s/prism.ini" % (workPath)
    if not (os.path.exists(configPath)):
    	shutil.copy2(os.path.join(run_files_path, "prism.ini"), configPath)
    
    Controller(jobId)
else:
	print "usage: python prism.py <pdbListPath> <templateListPath> <jobId>"
