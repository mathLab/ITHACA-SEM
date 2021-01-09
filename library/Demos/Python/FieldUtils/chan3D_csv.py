#-f -e -m interppoints:fromxml=chan3D.xml:fromfld=chan3D.fld:topts=chan3D_pts.csv out.csv
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

ProcessModule.Create("interppoints", field, fromxml="chan3D.xml", fromfld="chan3D.fld", topts="chan3D_pts.csv").Run()
OutputModule.Create("csv", field, "out.csv", test=True).Run()
