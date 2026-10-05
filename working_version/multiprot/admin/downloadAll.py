#!/usr/bin/python
import os
import MySQLdb as mdb
import datetime

i = datetime.datetime.now()
fileName = "/home/prism/projects/prism/results/%d_%d_%d.txt" % (i.year,i.month,i.day)
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

#clean previous file
os.system("rm /home/prism/projects/prism/results/*")

query = "select target1,target2,interface,energy,date_column from results order by date_column desc"

cur.execute(query)

rows = cur.fetchall()

filehnd = open(fileName,"w")
filehnd.write("PRISM Web Server Predictions:\n")
filehnd.write("Last Update: %d/%d/%d\n" % (i.year,i.month,i.day))
filehnd.write("File Format: Target1 Target2 Interface Energy Date_Added\n\n")
for row in rows:
	filehnd.write("%s\t%s\t%s\t%s\t%s\n" % (row[0],row[1],row[2],row[3],row[4]))

filehnd.close()
if con:
	con.close()
