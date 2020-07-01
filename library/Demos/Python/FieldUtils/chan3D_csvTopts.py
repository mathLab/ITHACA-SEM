# -f -e chan3D_pts.csv chan3D_pts.pts
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

InputModule.Create("csv", field, infile="chan3D_pts.csv", addfiles="csv:chan3D_pts.csv").Run()
OutputModule.Create("pts" , field, test=True, outfile="chan3D_pts.pts").Run()
