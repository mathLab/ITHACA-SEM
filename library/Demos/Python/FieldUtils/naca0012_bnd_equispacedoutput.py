# -f -e -m equispacedoutput naca0012_bnd.xml naca0012_b0.fld output.vtu
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

InputXml(field, "naca0012_bnd.xml").Run()
InputFld(field, "naca0012_b0.fld").Run()
ProcessEquiSpacedOutput(field).Run()
OutputVtk(field, "output.vtu").Run()

