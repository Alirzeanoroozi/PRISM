#!/usr/bin/env python
# Attila Gursoy, Oct 7 2016
import sys
import subprocess

#  convert torque jobid to prism jobid
#  usage: python queryJob.py  <torque jobid>

jobsFolder = "/home/prism/projects/prism/jobs/"
if (len(sys.argv) != 2):
   print "usage: python torque2prism.py  <torque jobid>"
   exit()


torquejob = sys.argv[1]
output = subprocess.check_output("qstat "+torquejob, shell=True)
jobsuffix=output.split()[15][8:]
output = subprocess.check_output("ls -d {}*{}".format(jobsFolder,jobsuffix), shell=True)
print output.split('/')[-1]
