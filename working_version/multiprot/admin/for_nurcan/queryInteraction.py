import shutil
import os
import MySQLdb as mdb

prismDir = "/home/prism/projects/prism/"
modelDir = "./models/"

db_f = open("../../config.inc", "r")
db_f.readline()
my_host = db_f.readline().split("'")[1]
my_user = db_f.readline().split("'")[1]
my_pass = db_f.readline().split("'")[1]
my_db = db_f.readline().split("'")[1]
db_f.close()
dbcon = mdb.connect(host=my_host, user=my_user, passwd=my_pass, db=my_db)



cur = dbcon.cursor()
cur.execute("SELECT VERSION()")
ver = cur.fetchone()

print "Database version : %s " % ver



def   searchInteraction(p1,p2,interface,energy,line):
    cur = dbcon.cursor()

    cur.execute("SELECT * FROM results where (BINARY target1=%s && BINARY target2=%s && BINARY interface=%s && BINARY energy=%s)",(p1,p2,interface,energy))
    rows = cur.fetchall()
    for row in rows: 
        #print row
        source_pdb = prismDir+row[4]
        source_txt = source_pdb[0:-13]+"intRes.txt"
        #print source_pdb, source_txt
        #print row[0], row[1], row[2], row[3], row[4].split("/")[-1]
        shutil.copy(source_pdb,modelDir)
        shutil.copy(source_txt,modelDir)
        print line+" "+source_pdb.split("/")[-1][0:-14]
     

with open("Mutated_PRISM_Interaction.txt") as f:
    for line in f:
        words = line.strip().split()
        if words[4] != "None":
           words[4] = words[4]+" "+words[5] 
        p1 = words[0]
        p2 = words[1]
        interface = words[2]
        energy = words[3]
        date = words[4]
         
        #searchInteraction(p1.lower(),p2.lower());
        searchInteraction(p1,p2, interface, energy,line.strip())



dbcon.close();
