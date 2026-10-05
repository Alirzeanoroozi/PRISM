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
	cur.execute("SELECT * FROM templates where template=%s",(temp[1].strip()))
	rows = cur.fetchall()
	if len(rows) != 1:
		print line				
	else:
		print line
if con:
	con.close()
fi.close()
