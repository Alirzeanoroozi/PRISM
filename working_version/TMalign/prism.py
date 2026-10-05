#!/usr/bin/env python2.7
import os
import runpy


root = os.path.abspath(os.path.dirname(__file__))
os.chdir(root)
runpy.run_path(os.path.join(root, "run_files", "prism.py"), run_name="__main__")
