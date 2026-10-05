#!/usr/bin/python
import os
import MySQLdb as mdb

jobsFolder = "/home/prism/projects/prism/jobs/"
fiberdockFolder = "/home/prism/projects/prism/fiberdock_output/"
pdbFolder = "/home/prism/projects/prism/pdb/"
#this scripts first checks jobs that are older than 7 days
(a,b) = os.popen4("find %s -type d -mtime +7" % (jobsFolder))
k = b.readlines()
a.close()
b.close()
#database connect
db_f = open("../config.inc", "r")
db_f.readline()
my_host = db_f.readline().split("'")[1]
my_user = db_f.readline().split("'")[1]
my_pass = db_f.readline().split("'")[1]
my_db = db_f.readline().split("'")[1]
db_f.close()
con = mdb.connect(host=my_host, user=my_user, passwd=my_pass, db=my_db)
cur = con.cursor()
for line in k:
	line = line.strip().split(jobsFolder)[1]
	#delete jobs from the jobs folder
	job = line[:26]
	#check if the job is still in the queue
	cur.execute("select job_id from ip_addr where job_id = %s",(job))
	row = cur.fetchall()
	if len(row) == 0:
		os.system("rm -rf %s%s" % (jobsFolder,job))
		#fetch files from database
		cur.execute("select structure from results where structure like %s",("%"+job+"%"))
		rows = cur.fetchall()
		if len(rows) == 0:
			os.system("rm -rf %s%s" % (fiberdockFolder,job)) #clean fiberdock_output folder
#clean pdb folder
#os.system("rm -rf %s*" % (pdbFolder))

if con:
	con.close()	
