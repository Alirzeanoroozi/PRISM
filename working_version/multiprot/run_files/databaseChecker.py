
#!/usr/bin/env python
#Written by Alper Baspinar
#creates mysql queries
import MySQLdb as mdb
class DatabaseChecker:
	global leftTarget,rightTarget,cur,con
	#constructor of the class requires leftTarget and rightTarget lists
	def __init__(self,leftTarget,rightTarget):
		db_f = open("../config.inc", "r")
		db_f.readline()
		my_host = db_f.readline().split("'")[1]
		my_user = db_f.readline().split("'")[1]
		my_pass = db_f.readline().split("'")[1]
		my_db = db_f.readline().split("'")[1]
		db_f.close()
		self.con = mdb.connect(host=my_host, user=my_user, passwd=my_pass, db=my_db)
		self.cur = self.con.cursor()

		self.leftTarget = leftTarget
		self.rightTarget = rightTarget
		
	
	def checker(self):
		preLeftTarget = []
		preRightTarget = []
		tableEntry = []
		for index in range(len(self.leftTarget)):
			pdb1 = self.leftTarget[index]
			pdb2 = self.rightTarget[index]
			#if len(self.templateList) == 1:
			#	self.cur.execute("SELECT * FROM results where BINARY target1=%s && BINARY target2=%s && BINARY interface=%s order by energy ASC",(pdb1,pdb2,self.templateList[0]))
			#else:
			self.cur.execute("SELECT * FROM results where (BINARY target1=%s && BINARY target2=%s) || (BINARY target1=%s && BINARY target2=%s) order by energy ASC",(pdb1,pdb2,pdb2,pdb1))
			rows = self.cur.fetchall()
			if len(rows) == 0:
				self.cur.execute("SELECT * FROM jobs where BINARY target1=%s && BINARY target2=%s",(pdb1,pdb2))
				if len(rows) == 0:
					self.cur.execute("SELECT * FROM jobs where BINARY target1=%s && BINARY target2=%s",(pdb2,pdb1))
					if len(rows) == 0:
						preLeftTarget.append(pdb1)
						preRightTarget.append(pdb2)
			else:
				tableEntry.append(rows)
				 		
	
		if self.con:
			self.con.close()
		
		return preLeftTarget,preRightTarget,tableEntry
