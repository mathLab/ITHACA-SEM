# -f -m vorticity -e chan3D.xml chan3D.fld chan3D_vort.fld
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

InputModule.Create("xml", field, "chan3D.xml").Run()
InputModule.Create("fld", field, "chan3D.fld").Run()
ProcessModule.Create("vorticity", field).Run()
OutputModule.Create("fld", field, "chan3D_vort.fld").Run()
