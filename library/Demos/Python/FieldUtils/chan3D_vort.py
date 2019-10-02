# -f -m vorticity -e chan3D.xml chan3D.fld chan3D_vort.fld
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

InputXml(field, "chan3D.xml").Run()
InputFld(field, "chan3D.fld").Run()
ProcessVorticity(field).Run()
OutputFld(field, "chan3D_vort.fld").Run()
