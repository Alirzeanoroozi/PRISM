#!/usr/bin/env python
#Written by Alper Baspinar
#creates mysql queries
import MySQLdb as mdb

#constructor of the class requires mode and the list to write to database
my_host = "localhost"
my_user = "prismUser"
my_pass = "prismPass"
my_db = "prismDatabase"
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
