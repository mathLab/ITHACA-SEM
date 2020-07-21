# -f -e -m fieldfromstring:fieldstr="p*x+y":fieldname="p" chan3D.xml chan3D.fld chan3D_modp.csv
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

InputModule.Create("xml", field, "chan3D.xml").Run()
InputModule.Create("fld", field, "chan3D.fld").Run()
ProcessModule.Create("fieldfromstring", field, fieldstr="p*x+y", fieldname="p").Run()
OutputModule.Create("pts", field, "chan3D_modp.csv").Run()
