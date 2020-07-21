# -f -e -n 10 chan3D.xml chan3D.fld chan3D.dat
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True, output_points=10)

print("hello pytest")
InputModule.Create("xml",  field, "chan3D.xml").Run()
InputModule.Create("fld",  field, "chan3D.fld").Run()
OutputModule.Create("dat", field, "chan3D.dat").Run()
