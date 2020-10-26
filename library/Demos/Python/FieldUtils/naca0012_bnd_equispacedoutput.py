# -f -e -m equispacedoutput naca0012_bnd.xml naca0012_b0.fld output.vtu
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

InputModule.Create("xml", field, infile={"xml":"naca0012_bnd.xml"}).Run()
InputModule.Create("fld", field, infile={"fld":"naca0012_b0.fld"}).Run()
ProcessModule.Create("equispacedoutput", field).Run()
OutputModule.Create("vtu", field, outfile="output.vtu").Run()
