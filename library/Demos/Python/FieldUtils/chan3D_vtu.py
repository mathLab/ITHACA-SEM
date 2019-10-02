# -f chan3D.xml chan3D.fld chan3D.vtu
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True)

InputXml(field, "chan3D.xml").Run()
InputFld(field, "chan3D.fld").Run()
OutputVtk(field, "chan3D.vtu").Run()
