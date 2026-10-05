#!/usr/bin/env python
#Written by Alper Baspinar
#send mail
import os,sys
import MySQLdb as mdb
import smtplib
jobId = sys.argv[1]
my_host = "localhost"
my_user = "prismUser"
my_pass = "prismPass"
my_db = "prismDatabase"
con = mdb.connect(host=my_host, user=my_user, passwd=my_pass, db=my_db)
cur = con.cursor()
cur.execute("select * from job where jobid=%s",(jobId,))
rows = cur.fetchall()
if len(rows) != 0 and rows[0][3] != "":
	fromaddr = 'prism@ku.edu.tr'
	toaddrs  = rows[0][3]
	msg = "\r\n".join([
	  "From: no-reply@ku.edu.tr",
	  "Reply-To: no-reply@ku.edu.tr",
	  "To: %s" % (rows[0][3]),
	  "Subject: Prism-Webserver Run: %s is submitted."% (jobId),
	  "",
	  "Your job is submitted successfully. You will receive another email when your Prism run is completed. You can bookmark the link given below. Thanks for using Prism.\n\nhttp://cosbi.ku.edu.tr/prism/result.php?jobId=%s\n" % (jobId)])

	username = 'prism@ku.edu.tr'
	password = 'fehjnmnpcmbnealo'
	server = smtplib.SMTP('smtp.gmail.com:587')
	server.ehlo()
	server.starttls()
	server.login(username,password)
	server.sendmail(fromaddr, toaddrs, msg)
	server.quit()
if con:
	con.close()
