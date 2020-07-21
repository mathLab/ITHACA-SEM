# -f -e chan3D_pts.csv chan3D_pts.pts
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

InputModule.Create("csv",   field, "chan3D_pts.csv").Run()
OutputModule.Create("pts" , field, "chan3D_pts.pts", test=True).Run()
