#!/usr/bin/env python
# Attila Gursoy, Oct 7, 2016
import os
import sys
import MySQLdb as mdb

#  deletes a job whose jobid is provided in the command line
#  usage: deleteJob.py  <jobid>

jobsFolder = "/home/prism/projects/prism/jobs/"
fiberdockFolder = "/home/prism/projects/prism/fiberdock_output/"
pdbFolder = "/home/prism/projects/prism/pdb/"

if (len(sys.argv) != 2):
   print "usage: python deleteJob.py  <jobid>"
   exit()

job = sys.argv[1]
print "Removing job " + job

ip = "NULL"
jobid = "NULL"
target1 = "NULL"
target2 = "NULL"
email = "NULL"
jobdate = "NULL"

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

#delete jobs from the jobs folder that are running but stucked

#check if the job is still in the queue
cur.execute("select * from ip_addr where job_id = %s",(job))
row = cur.fetchall()

if len(row) == 1:
   (ip,jobid) = row[0]
   
   print "IP " + ip + " Jobid " + jobid
   print "remove job from ip_addr table"
   print "IP " + row[0][0] + " Jobid" + row[0][1]
   cur.execute("delete from ip_addr where job_id = %s",(job))
   

   print "remove job from job table"

   cur.execute("select * from job where jobid = %s",(job))
   jobrow = cur.fetchall()
   if len(jobrow) == 1:
        (tmp,target1,target2,email,date) = jobrow[0]
        print  "remove ", target1,target2,email,date
        cur.execute("delete from job where jobid = %s",(job))
   
   cur.execute("select structure from results where structure like %s",("%"+job+"%"))
   fiberdock_rows = cur.fetchall()
   if  len(fiberdock_rows) == 0:
       # no fiberdock results entered into database, we can remove fiberdock directory
       print "no results in database, remove partial results from fiberdock_output directory if there is"
       os.system("rm -rf %s%s" % (fiberdockFolder,job)) #clean fiberdock_output folder
   else:
       print "results table has data, please check fiberdock_output for consistency"

   # remove job folder
   print "remove job folder"
   os.system("rm -rf %s%s" % (jobsFolder,job))

if con:
        con.commit()
        con.close()

print jobid, ip, target1, target2, email, date
