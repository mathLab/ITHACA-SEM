# -f -e -m fieldfromstring:fieldstr="p*x+y":fieldname="p" chan3D.xml chan3D.fld chan3D_modp.csv
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

InputXml(field, "chan3D.xml").Run()
InputFld(field, "chan3D.fld").Run()
ProcessFieldFromString(field, fieldstr="p*x+y", fieldname="p").Run()
OutputPts(field, "chan3D_modp.csv").Run()
