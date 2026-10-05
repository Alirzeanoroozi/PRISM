#!/usr/bin/env python
# Attila Gursoy, Oct 7 2016
import os
import sys
import MySQLdb as mdb

#  query a job whose jobid is provided in the command line
#  usage: python queryJob.py  <jobid>

jobsFolder = "/home/prism/projects/prism/jobs/"
fiberdockFolder = "/home/prism/projects/prism/fiberdock_output/"
pdbFolder = "/home/prism/projects/prism/pdb/"

if (len(sys.argv) != 2):
   print "usage: queryJob.py  <jobid>"
   exit()

job = sys.argv[1]

ip = "NULL"
jobid = "NULL"
target1 = "NULL"
target2 = "NULL"
email = "NULL"
jobdate = "NULL"

print
print "Query job queue for",job
os.system("qstat | grep %s"%(job[18:]))
print

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

#check if the job is still in the queue
cur.execute("select * from ip_addr where job_id = %s",(job))
row = cur.fetchall()

if len(row) == 1:
   (ip,jobid) = row[0]
   
   print "table ip_addr: IP " + ip + " Jobid " + jobid

   cur.execute("select * from job where jobid = %s",(job))
   jobrow = cur.fetchall()
   if len(jobrow) == 1:
        (tmp,target1,target2,email,date) = jobrow[0]
        print "job table ",tmp,target1,target2, email,date
  
   print 
   print "fiberdock_output folder:"
   print "ls -al {}{}".format(fiberdockFolder,job)
   os.system("ls -al %s%s" % (fiberdockFolder,job))
   cur.execute("select structure from results where structure like %s",("%"+job+"%"))
   fiberdock_rows = cur.fetchall()
   if  len(fiberdock_rows) == 0:
       print "no results in database"
   else:
       print "results table:"
       print fibderdock_rows
 
   print 
   print "job folder:"
   print "ls -al {}{}".format(jobsFolder,job)
   os.system("ls -al {}{}".format(jobsFolder,job))

if con:
        con.close()
