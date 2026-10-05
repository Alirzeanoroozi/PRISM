#!/usr/bin/env python
#Written by Alper Baspinar
#send mail
import os
import MySQLdb as mdb
import smtplib
class MailSender:
	#constructor of the class requires jobId
	def __init__(self,jobId):
		#workPath = "../jobs/"+jobId
		my_host = "localhost"
		my_user = "prismUser"
		my_pass = "prismPass"
		my_db = "prismDatabase"
		self.con = mdb.connect(host=my_host, user=my_user, passwd=my_pass, db=my_db)
		self.cur = self.con.cursor()
		self.cur.execute("select * from job where jobid=%s",(jobId,))
		rows = self.cur.fetchall()
		if len(rows) != 0 and rows[0][3] != "":
			fromaddr = 'prism@ku.edu.tr'
			toaddrs  = rows[0][3]
			msg = "\r\n".join([
			  "From: no-reply@ku.edu.tr",
			  "Reply-To: no-reply@ku.edu.tr",
			  "To: %s" % (rows[0][3]),
			  "Subject: Prism-Webserver Run: %s is completed."% (jobId),
			  "",
			  "Your Prism run is completed. You can reach the results using the link given below. The link will be deleted after 7 days. Thanks for using Prism.\n\nhttp://cosbi.ku.edu.tr/prism/result.php?jobId=%s\n" % (jobId)])

			username = 'prism@ku.edu.tr'
			password = 'fehjnmnpcmbnealo'
			server = smtplib.SMTP('smtp.gmail.com:587')
			server.ehlo()
			server.starttls()
			server.login(username,password)
			server.sendmail(fromaddr, toaddrs, msg)
			server.quit()
		if self.con:
                        self.con.close()
