# -f chan3D.xml chan3D.fld chan3D.vtu
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True)

InputModule.Create("xml", field, infile={"xml":"chan3D.xml"}).Run()
InputModule.Create("fld", field, infile={"fld":"chan3D.fld"}).Run()
OutputModule.Create("vtu", field, outfile="chan3D.vtu").Run()
