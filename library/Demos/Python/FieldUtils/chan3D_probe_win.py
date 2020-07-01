# -f -m interppoints:fromxml=chan3D.xml:fromfld=chan3D.fld chan3D_probe.pts chan3D_probe.dat
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, orceoutput=True)

ProcessModule.Create("interppoints", field, fromxml="chan3D.xml", fromfld="chan3D.fld", topts="chan3D_probe.pts").Run()
OutputModule.Create("dat", field, outfile="chan3D_probe.dat").Run()
