#!/usr/bin/env python
#Written by Alper Baspinar
#creates mysql queries
import MySQLdb as mdb

class MysqlWriter:
	#constructor of the class requires mode and the list to write to database
	def __init__(self,mode,aList):
		db_f = open("../config.inc", "r")
		db_f.readline()
		my_host = db_f.readline().split("'")[1]
		my_user = db_f.readline().split("'")[1]
		my_pass = db_f.readline().split("'")[1]
		my_db = db_f.readline().split("'")[1]
		db_f.close()
		con = mdb.connect(host=my_host, user=my_user, passwd=my_pass, db=my_db)
		cur = con.cursor()
		size = len(aList)
		if mode == 0 and size !=0:
			for i in aList:
				if i != "pdb1" and i != "pdb2":
					try:
						cur.execute("INSERT INTO targets (target) VALUES (%s)", (i))
						con.commit()
					except:
						con.rollback()
		
		elif mode == 1 and size!=0:
			leftTarget = aList[0]
			rightTarget = aList[1]
			for i in range(len(leftTarget)):
				left = leftTarget[i]
				right = rightTarget[i]
				if left != "pdb1" and left != "pdb2" and right != "pdb1" and right != "pdb2":
					try:
						cur.execute("INSERT INTO jobs (target1,target2) VALUES (%s,%s)", (left,right))
						con.commit()
					except:
						con.rollback()
		
		elif mode == 2 and size!=0:
			for i in aList:
				temp = i.split("_")
				if temp[2] != "pdb1" and temp[2] != "pdb1":
					try:
						cur.execute("INSERT INTO passed (target,interface,chain) VALUES (%s,%s,%s)", (temp[2],temp[0],temp[1]))
						con.commit()
					except:
						con.rollback()
		
		elif mode == 3 and size!=0:
			for i in aList:
				e = i[0]
				structure = i[1]
				temp = e.split()
				energy = temp[2]
				a = temp[0].split("_")
				interface = a[0]
				target1 = a[2].split(".")[0]
				a = temp[1].split("_")
				target2 = a[2].split(".")[0]
				if target1 != "pdb1" and target1 != "pdb2" and target2 != "pdb1" and target2 != "pdb2":
					try:
						cur.execute("INSERT INTO results (target1,target2,interface,energy,structure) VALUES (%s,%s,%s,%s,%s)", (target1,target2,interface,energy,structure))
						con.commit()
					except:
						con.rollback()
		elif mode == 5:
			try:
				cur.execute("Delete from ip_addr where job_id=%s",(aList))
				con.commit()
			except:
				con.rollback()
					
		if con:
			con.close()

