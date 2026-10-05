import os,sys
path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'run_files'))
if not path in sys.path:
    sys.path.insert(1, path)
from mainController import Controller

if len(sys.argv) == 4:
    pdbListPath = sys.argv[1]
    templateListPath = sys.argv[2]
    jobId = sys.argv[3]
    currentPath = os.getcwd()
    #creates a folder for the workPath
    workPath = "jobs/%s" % jobId
    if os.path.exists(workPath) == False:
	os.mkdir(workPath,0777)
    #creates a folder for the lists
    listPath = workPath+"/lists"
    if os.path.exists(listPath) == False:
        os.mkdir(listPath,0777)
    #moves pdbList to workFolder
    pdbPath = listPath+"/pair_list"
    if not (os.path.exists(pdbPath)):
        os.system("cp %s %s" % (pdbListPath,pdbPath))
    #moves templateList to workFolder
    templatePath = listPath+"/template_list"
    if not (os.path.exists(templatePath)):
        os.system("cp %s %s" % (templateListPath,templatePath))
    
    configPath = "%s/prism.ini" % (workPath)
    if not (os.path.exists(configPath)):
    	os.system("cp %s %s" % ("prism.ini",configPath))
    
    Controller(jobId)
else:
	print "usage: python prism.py <pdbListPath> <templateListPath> <jobId>"
	
del path
