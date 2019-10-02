# -f -e -m interppoints:fromxml=interpfieldline.xml:fromfld=interpfieldline.fld:line=10,0.024,0.0,0.16,0.0 interpfieldline.dat
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

ProcessInterpPoints(field, fromxml="interpfieldline.xml", fromfld="interpfieldline.fld", line="10,0.024,0.0,0.16,0.0").Run()
OutputTecplot(field, "interpfieldline.dat").Run()
