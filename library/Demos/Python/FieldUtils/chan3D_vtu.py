# -f chan3D.xml chan3D.fld chan3D.vtu
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True)

InputModule.Create("xml",  field, "chan3D.xml").Run()
InputModule.Create("fld",  field, "chan3D.fld").Run()
OutputModule.Create("vtu", field, "chan3D.vtu").Run()
