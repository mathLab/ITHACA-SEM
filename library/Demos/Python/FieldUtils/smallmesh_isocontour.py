# -f -e -n6 -m isocontour:fieldstr="p":fieldvalue=0.1:globalcondense:smooth smallmesh.xml smallmesh.fld iso.dat
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True, output_points=6)

print("input modules")
InputModule.Create("xml", field, "smallmesh.xml").Run()
InputModule.Create("fld", field, "smallmesh.fld").Run()
print("process modules")
ProcessModule.Create("equispacedoutput", field).Run()
ProcessModule.Create("isocontour", field, fieldstr="p", fieldvalue="0.1", globalcondense="1", smooth="1").Run()
print("output modules")
OutputModule.Create("dat", field, "iso.dat").Run()
