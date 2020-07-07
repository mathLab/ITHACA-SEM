# -f -m vorticity -e chan3D.xml chan3D.fld chan3D_vort.fld
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

InputModule.Create("xml", field, infile={"xml":"chan3D.xml"}).Run()
InputModule.Create("fld", field, infile={"fld":"chan3D.fld"}).Run()
ProcessModule.Create("vorticity", field).Run()
OutputModule.Create("fld", field, outfile="chan3D_vort.fld").Run()
