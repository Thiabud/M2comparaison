# -*- coding: utf-8 -*-
#!/usr/bin/python3
import os
import sys
import subprocess
import shutil
import argparse
import csv
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from collections import Counter
#----------------------------------#
#    PARAMETRES
#----------------------------------#

#GESTION DES ARGUMENTS 
parser = argparse.ArgumentParser(description='')
parser.add_argument("-r", "--ref", dest="ref" , type=str, help="directory with targert region of ref ") 
parser.add_argument("-n", "--new", dest="new" , type=str, help="directory with target region of new") 

args = parser.parse_args()
new = args.new
ref = args.ref

countnumberoftargetregioninref = "ls " + ref + " | wc -l"
countnumberoftargetregioninnew = "ls " + new + " | wc -l"
targetinref = subprocess.check_output(countnumberoftargetregioninref, shell=True)
targetinnew = subprocess.check_output(countnumberoftargetregioninnew, shell=True)
targetinref=int(targetinref)
targetinnew=int(targetinnew)

print("targetinref:", targetinref)
print("targetinnew:", targetinnew)

numberofnuclinref="cat " + ref + "* | wc -c"
numberofnuclinref = subprocess.check_output(numberofnuclinref, shell=True)
numberofnuclinref=int(numberofnuclinref)
print("numberofnuclinref",numberofnuclinref)

numberofnuclinnew="cat " + new + "* | wc -c"
numberofnuclinnew = subprocess.check_output(numberofnuclinnew, shell=True)
numberofnuclinnew=int(numberofnuclinnew)
print("numberofnuclinnew",numberofnuclinnew)