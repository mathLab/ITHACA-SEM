# -f -e -m isocontour:fieldstr="u+v":fieldvalue=0.5:fieldname="UplusV":globalcondense:smooth chan3D.xml chan3D.fld isocontour.dat
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

InputModule.Create("xml", field, "chan3D.xml").Run()
InputModule.Create("fld", field, "chan3D.fld").Run()
ProcessModule.Create("equispacedoutput", field).Run()
ProcessModule.Create("isocontour", field, fieldstr="u+v", fieldvalue="0.5", fieldname="UplusV", globalcondense="1", smooth="1").Run()
OutputModule.Create("dat", field, "isocontour.dat").Run()
