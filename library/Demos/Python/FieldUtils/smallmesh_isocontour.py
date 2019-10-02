# -f -e -n6 -m isocontour:fieldstr="p":fieldvalue=0.1:globalcondense:smooth smallmesh.xml smallmesh.fld iso.dat
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True, output_points=6)

InputXml(field, "smallmesh.xml").Run()
InputFld(field, "smallmesh.fld").Run()
ProcessEquiSpacedOutput(field).Run()
ProcessIsoContour(field, fieldstr="p", fieldvalue="0.1", globalcondense=True, smooth=True).Run()
OutputTecplot(field, "iso.dat").Run()
