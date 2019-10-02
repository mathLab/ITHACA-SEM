# -f -m equispacedoutput chan3D.xml chan3D.fld equispacedoutput.dat
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True)

InputXml(field, "chan3D.xml").Run()
InputFld(field, "chan3D.fld").Run()
ProcessEquiSpacedOutput(field).Run()
OutputTecplot(field, "equispacedoutput.dat").Run()
