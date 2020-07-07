# -f -m meanmode -e chan3DH1D.xml chan3DH1D.fld chan3DH1D_mean.fld
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

InputModule.Create("xml", field, infile={"xml":"chan3DH1D.xml"}).Run()
InputModule.Create("fld", field, infile={"fld":"chan3DH1D.fld"}).Run()
ProcessModule.Create("meanmode", field).Run()
OutputModule.Create("fld", field, outfile="chan3DH1D_mean.fld").Run()
