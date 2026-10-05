#!/usr/bin/python
import os
import MySQLdb as mdb
import smtplib

def sendMail(jobId):
	fromaddr = 'prism@ku.edu.tr'
        toaddrs  = 'prism@ku.edu.tr'
        msg = "\r\n".join([
          "From: no-reply@ku.edu.tr",
          "Reply-To: no-reply@ku.edu.tr",
          "To: %s" % ("prism@ku.edu.tr"),
          "Subject: Job %s is re-queued."% (jobId),
          "",
          "Hi! Dear Master, \n\nThere was a problem with the job %s. The reason might be the server is shut down or there might be a bug with prism webserver. This message is sent as a result of the regular check from crontab. Please keep an eye on the progress of this job. Thank you. \n\n-alper." % (jobId)])

        username = 'prism@ku.edu.tr'
        password = 'fehjnmnpcmbnealo'
        server = smtplib.SMTP('smtp.gmail.com:587')
        server.ehlo()
        server.starttls()
        server.login(username,password)
        server.sendmail(fromaddr, toaddrs, msg)
        server.quit()

#jobs folder
jobs = "/home/prism/projects/prism/jobs"
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
cur.execute("select job_id from ip_addr")
rows = cur.fetchall()
if len(rows) != 0:
	for row in rows:
		(a,b) = os.popen4("qstat -f | grep \"Job_Name = %s\"" % (row[0]))
		k = b.readlines()
		a.close()
		b.close()
		if len(k) == 0:
			os.system("/usr/bin/touch %s/%s" % (jobs,row[0]))
			os.system("/usr/bin/qsub -N %s -v job_id=%s /home/prism/projects/prism/prism.sh" % (row[0],row[0]))
			sendMail(row[0])
			
if con:
	con.close()

