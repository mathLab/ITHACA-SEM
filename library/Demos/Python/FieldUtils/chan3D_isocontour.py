# -f -e -m isocontour:fieldstr="u+v":fieldvalue=0.5:fieldname="UplusV":globalcondense:smooth chan3D.xml chan3D.fld isocontour.dat
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

InputXml(field, "chan3D.xml").Run()
InputFld(field, "chan3D.fld").Run()
ProcessEquiSpacedOutput(field).Run()
ProcessIsoContour(field, fieldstr="u+v", fieldvalue="0.5", fieldname="UplusV", globalcondense=True, smooth=True).Run()
OutputTecplot(field, "isocontour.dat").Run()
