# -f -e -m interppoints:fromxml=interpfieldline.xml:fromfld=interpfieldline.fld:line=10,0.024,0.0,0.16,0.0 interpfieldline.dat
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

ProcessModule.Create("interppoints", field, fromxml="interpfieldline.xml", fromfld="interpfieldline.fld", line="10,0.024,0.0,0.16,0.0").Run()
OutputModule.Create("dat", field, outfile="interpfieldline.dat").Run()
