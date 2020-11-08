# -f -e chan3D_pts.pts chan3D_pts.csv
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

InputModule.Create("pts", field, "chan3D_pts.pts").Run()
OutputModule.Create("csv", field, "chan3D_pts.csv").Run()
