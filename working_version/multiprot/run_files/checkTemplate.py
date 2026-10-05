#!/usr/bin/env python
#written by Alper Baspinar
#checks the template and creates if necessary
from templateGenerator import TemplateGenerator
class TemplateChecker:
	#constructor of the class requires templateList
	global workPath,templateList
	def __init__(self,workPath,templateList):
		self.templateList = templateList
		self.workPath = workPath
	def checker(self):
		#read default templateList
		filehnd = open("../template_default","r")
		allTemplates = []
		for line in filehnd.readlines():
			line = line.strip()[:6]
			allTemplates.append(line)
		noTemplateList = []
		yesTemplateList = []
		for temp in self.templateList:
			if not (temp in allTemplates):
				noTemplateList.append(temp)
			else:
				yesTemplateList.append(temp)
		if len(noTemplateList) == 0:
			size = len(yesTemplateList)
			if size == len(allTemplates):
				return [1,self.templateList]
			else:
				return [2,yesTemplateList]
		else:
			#try to generate template files
			gen = TemplateGenerator(self.workPath,noTemplateList).generator()
			if gen[0] == 0:
				size = len(yesTemplateList)
				if size == 0:
					return [0,[]]
				else:
					#check if yesTemplateList == default template
					if size == len(allTemplates):
						return [1,yesTemplateList]
					else:
						return [2,yesTemplateList]
			else:
				return [2,yesTemplateList+gen[1]]
		filehnd.close()		
