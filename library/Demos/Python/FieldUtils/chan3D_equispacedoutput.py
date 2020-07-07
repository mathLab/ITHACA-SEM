# -f -m equispacedoutput -e chan3D.xml chan3D.fld equispacedoutput.dat
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

InputModule.Create("xml", field, infile={"xml":"chan3D.xml"}).Run()
InputModule.Create("fld", field, infile={"fld":"chan3D.fld"}).Run()
ProcessModule.Create("equispacedoutput", field).Run()
OutputModule.Create("dat", field, outfile="equispacedoutput.dat").Run()
