# -f -e -n6 -m isocontour:fieldstr="p":fieldvalue=0.1:globalcondense:smooth smallmesh.xml smallmesh.fld iso.dat
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True, output_points=6)

InputModule.Create("xml", field, infile="smallmesh.xml", addfiles="xml:smallmesh.xml").Run()
InputModule.Create("fld", field, infile="smallmesh.fld", addfiles="fld:smallmesh.fld").Run()
ProcessModule.Create("equispacedoutput", field).Run()
ProcessModule.Create("isocontour", field, fieldstr="p", fieldvalue="0.1", globalcondense=True, smooth=True).Run()
OutputModule.Create("dat", field, outfile="iso.dat").Run()
