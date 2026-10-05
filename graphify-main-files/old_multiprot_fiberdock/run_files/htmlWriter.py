#!/usr/bin/env python
#Written by Alper Baspinar
#creates html file for the job
class HtmlWriter:
    #constructor of the class requires filePath, mode and list to run
    def __init__(self,filePath,listPath,mode,aList):
	self.listPath = listPath
	listItem = ""
	size = len(aList)
	for i in aList:
		listItem += str(i) + " "
	if mode != 12 and mode != -1:
		filehnd = open(filePath,"a")
	else:
		filehnd = open(filePath,"w")
	
	if mode == 0:
		filehnd.write("Pdb download stage started...<br>")
	elif mode == -1:
		filehnd.write("<head><meta http-equiv=refresh content=60 > </head>");
                filehnd.write("<p><strong> The page will automatically refresh every minute. </strong></p>");
                filehnd.write("<p><strong> Job Status: </strong>Running</p>");
	elif mode == 1:
		if size == 0:
			filehnd.write("Check pdb inputs again!!!<br>")
		else:
			filehnd.write("%sdownloaded...<br>" % listItem)
		filehnd.write("<script>$progressUnit += 10;</script>")
	elif mode == 2:
		filehnd.write("Pre-process stage started...<br>") 
	elif mode == 3:
		if size != 0: 
			filehnd.write("Pre-process for %scompleted...<br>" % listItem)
		filehnd.write("<script>$progressUnit += 10;</script>")
	elif mode == 4:
		if size != 0:
			filehnd.write("Surface extraction stage started...<br>")
	elif mode == 5:
		if size != 0:
			filehnd.write("Surface extraction for %scompleted...<br>" % listItem)
		filehnd.write("<script>$progressUnit += 10;</script>")
	elif mode == 6:
		if size != 0:
			filehnd.write("Structural alignment stage started...<br>")
	elif mode == 7:
		if size != 0:
			filehnd.write("Structural alignment for %scompleted...<br>" % listItem)
		filehnd.write("<script>$progressUnit += 30;</script>")
	elif mode == 8:
		if size != 0:
			filehnd.write("Transformation filtering stage started...<br>")
	elif mode == 9:
		if size != 0:
			filehnd.write("Transformation filtering for %scompleted...<br>" % listItem)
		filehnd.write("<script>$progressUnit += 10;</script>")
	elif mode == 10:
		if size != 0:
			filehnd.write("Flexible refinement stage started...<br>")
	elif mode == 11:
		if size != 0:
			filehnd.write("Flexible refinement for %scompleted...<br>" % listItem)
		filehnd.write("<script>$progressUnit += 30;</script>")
	elif mode == 12:
		size = len(aList[1])
		pdbLink = "http://www.rcsb.org/pdb/explore/explore.do?structureId="
                interfaceLink = "http://prism.ccbb.ku.edu.tr/hotregion/hotregionCrossLinkForSimpleRun.php?pdbName="
		if size != 0:
			a1 = self.infoWrite()
			for line in a1:
				filehnd.write(line)
			filehnd.write("<h4>Results:</h4>")
			filehnd.write("<table id=\"myTable\" class=\"table table-hover table-striped\"><thead><tr><th>Target1</th><th>Target2</th><th>Interface</th><th>Energy</th><th>Structure</th></tr></thead><tbody>")
			tableEntry = aList[0]
			aList = aList[1]
			energyList = []
			for i in aList:
				e = i[0]
				temp = e.split()
				energy = temp[2]
				try:
					energy = float(energy)
				except:
					energy = 1000000
				energyList.append(energy)
			energyList = sorted(range(len(energyList)),key=lambda k: energyList[k])
			
			for index in energyList:
				i = aList[index]
				e = i[0]
				structure = i[1]
				temp = e.split()
				energy = temp[2]
				try:
					energy = float(energy)
				except:
					energy = 1000000
				a = temp[0].split("_")
				interface = a[0]
				target1 = a[2].split(".")[0]
				a = temp[1].split("_")
				target2 = a[2].split(".")[0]
				llink = "#"
                                rlink = "#"
                                if target1[:4] == "pdb1" or target1[:4] == "pdb2":
                                	llink = "<a href=\"%s\" target=\"_blank\">%s</a>" % ("#",target1)
                                else:
                                        llink = "<a href=\"%s%s\" target=\"_blank\">%s</a>" % (pdbLink,target1[:4],target1)
                                if target2[:4] == "pdb1" or target2[:4] == "pdb2":
                                        rlink = "<a href=\"%s\" target=\"_blank\">%s</a>" % ("#",target2)
                                else:
                                       	rlink = "<a href=\"%s%s\" target=\"_blank\">%s</a>" % (pdbLink,target2[:4],target2)
				filehnd.write("<tr><td>%s</td><td>%s</td><td><a href=\"%s%s&chain_1=%s&chain_2=%s\" target=\"_blank\">%s</a></td><td>%s</td><td><a href=\"#myModal\" role=\"button\" class=\"btn\" data-toggle=\"modal\" onclick=\"setLink('%s','%s')\">View</a></td></tr>" % (llink,rlink,interfaceLink,interface[:4],interface[4],interface[5],interface,energy,structure,structure.split(".fiberdock.pdb")[0]+".intRes.txt"))
			#data from database
			for entry in tableEntry:                                        
				for e in entry:
					target1 = e[0]
                                        target2 = e[1]
                                        interface = e[2]
                                        energy = e[3]
                                        structure = e[4]
					llink = "#"
                                        rlink = "#"
                                        if target1[:4] == "pdb1" or target1[:4] == "pdb2":
						llink = "<a href=\"%s\" target=\"_blank\">%s</a>" % ("#",target1)
                                        else:
                                                llink = "<a href=\"%s%s\" target=\"_blank\">%s</a>" % (pdbLink,target1[:4],target1)
                                        if target2[:4] == "pdb1" or target2[:4] == "pdb2":
                                                rlink = "<a href=\"%s\" target=\"_blank\">%s</a>" % ("#",target2)
                                        else:
                                                rlink = "<a href=\"%s%s\" target=\"_blank\">%s</a>" % (pdbLink,target2[:4],target2)
                                	filehnd.write("<tr><td>%s</td><td>%s</td><td><a href=\"%s%s&chain_1=%s&chain_2=%s\" target=\"_blank\">%s</a></td><td>%s</td><td><a href=\"#myModal\" role=\"button\" class=\"btn\" data-toggle=\"modal\" onclick=\"setLink('%s','%s')\">View</a></td></tr>" % (llink,rlink,interfaceLink,interface[:4],interface[4],interface[5],interface,energy,structure,structure.split(".fiberdock.pdb")[0]+".intRes.txt"))
			filehnd.write("</tbody></table>")
		else:
			tableEntry = aList[0]
			if len(tableEntry) == 0:
				a1 = self.infoWrite()
                        	for line in a1:
                                	filehnd.write(line)
				filehnd.write("<h4>No Results Found</h4><br>")
			else:
				a1 = self.infoWrite()
                        	for line in a1:
                                	filehnd.write(line)
				filehnd.write("<h4>Results:</h4>")
                        	filehnd.write("<table id=\"myTable\" class=\"table table-hover table-striped\"><thead><tr><th>Target1</th><th>Target2</th><th>Interface</th><th>Energy</th><th>Structure</th></tr></thead><tbody>")
				for entry in tableEntry:
					for e in entry:
						target1 = e[0]
						target2 = e[1]
						interface = e[2]
						energy = e[3]
						structure = e[4]
						llink = "#"
						rlink = "#"
						if target1[:4] == "pdb1" or target1[:4] == "pdb2":
							llink = "<a href=\"%s\" target=\"_blank\">%s</a>" % ("#",target1)
						else:
							llink = "<a href=\"%s%s\" target=\"_blank\">%s</a>" % (pdbLink,target1[:4],target1)
						if target2[:4] == "pdb1" or target2[:4] == "pdb2":
							rlink = "<a href=\"%s\" target=\"_blank\">%s</a>" % ("#",target2)
						else:
							rlink = "<a href=\"%s%s\" target=\"_blank\">%s</a>" % (pdbLink,target2[:4],target2)
                                		filehnd.write("<tr><td>%s</td><td>%s</td><td><a href=\"%s%s&chain_1=%s&chain_2=%s\" target=\"_blank\">%s</a></td><td>%s</td><td><a href=\"#myModal\" role=\"button\" class=\"btn\" data-toggle=\"modal\" onclick=\"setLink('%s','%s')\">View</a></td></tr>" % (llink,rlink,interfaceLink,interface[:4],interface[4],interface[5],interface,energy,structure,structure.split(".fiberdock.pdb")[0]+".intRes.txt"))
				filehnd.write("</tbody></table>")
		
		filehnd.write("<script>$progressUnit = 100;</script>")
	else:
		filehnd.write("End...<br>")
	filehnd.close()
	
    def infoWrite(self):		
	a1 = ["<h4>Targets:</h4>"]
	filehnd = open("%s/pair_list" % (self.listPath))
	for line in filehnd.readlines():
		a1.append("<p>%s</p>" % (line))
	filehnd.close()
	a1.append("<h4>Template:</h4>")
	filehnd = open("%s/template_list" % (self.listPath))
	k = filehnd.readlines()
	if len(k) == 1:
		a1.append("<p>%s</p>" % (k[0]))
	else:
		a1.append("<p>%s</p>" % ("Default"))
	a1.append("<br>")
	filehnd.close()
        return a1			
			
		
