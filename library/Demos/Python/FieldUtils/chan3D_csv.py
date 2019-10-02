#-f -e -m interppoints:fromxml=chan3D.xml:fromfld=chan3D.fld:topts=chan3D_pts.csv out.csv
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

ProcessInterpPoints(field, fromxml="chan3D.xml", fromfld="chan3D.fld", topts="chan3D_pts.csv").Run()
OutputPts(field, "out.csv").Run()
