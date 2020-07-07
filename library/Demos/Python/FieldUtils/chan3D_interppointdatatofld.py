# -f -e -m interppointdatatofld:frompts=chan3D_pts.pts chan3D.xml chan3D_pts.fld
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

InputModule.Create("xml", field, infile={"xml":"chan3D.xml"}).Run()
ProcessModule.Create("interppointdatatofld", field, frompts="chan3D_pts.pts").Run()
OutputModule.Create("fld", field, outfile="chan3D_pts.fld").Run()
