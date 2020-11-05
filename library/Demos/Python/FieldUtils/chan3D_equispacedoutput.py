# -f -m equispacedoutput -e chan3D.xml chan3D.fld equispacedoutput.dat
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

InputModule.Create("xml",  field, "chan3D.xml").Run()
InputModule.Create("fld",  field, "chan3D.fld").Run()
ProcessModule.Create("equispacedoutput", field).Run()
OutputModule.Create("dat", field, "equispacedoutput.dat").Run()
