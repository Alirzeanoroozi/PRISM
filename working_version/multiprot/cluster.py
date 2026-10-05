#!/usr/bin/env python
#Written by Alper Baspinar
#creates mysql queries
import MySQLdb as mdb

#constructor of the class requires mode and the list to write to database
db_f = open("./config.inc", "r")
db_f.readline()
my_host = db_f.readline().split("'")[1]
my_user = db_f.readline().split("'")[1]
my_pass = db_f.readline().split("'")[1]
my_db = db_f.readline().split("'")[1]
db_f.close()
con = mdb.connect(host=my_host, user=my_user, passwd=my_pass, db=my_db)
cur = con.cursor()
fi = open("cluster_list")
for line in fi.readlines():
	line = line.strip()
	temp = line.split()
	try:
		cur.execute("INSERT INTO templates (cluster,template,repr) VALUES (%s,%s,%s)", (temp[0].strip(),temp[1].strip(),temp[2].strip()))
		con.commit()
	except:
		print "insert error for %s %s %s"  % (temp[0],temp[1],temp[2])
		con.rollback()
		
					
if con:
	con.close()
fi.close()
