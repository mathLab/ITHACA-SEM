# -f -e -m interppoints:fromxml=chan3D.xml:fromfld=chan3D.fld:topts=chan3D_probe.pts chan3D_probe.dat
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

ProcessInterpPoints(field, fromxml="chan3D.xml", fromfld="chan3D.fld", topts="chan3D_probe.pts").Run()
OutputTecplot(field, "chan3D_probe.dat").Run()
