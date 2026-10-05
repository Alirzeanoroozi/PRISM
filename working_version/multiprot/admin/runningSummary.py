import sys
import glob
lines = sys.stdin.readlines()

def jobSummary(line):
   job = line.split()[1]
   #print(job, job[8:])
   files = glob.glob("../jobs/*."+job[8:])
   print files[0][3:]


for line in lines:
    if "prism" in line:
       jobSummary(line)
