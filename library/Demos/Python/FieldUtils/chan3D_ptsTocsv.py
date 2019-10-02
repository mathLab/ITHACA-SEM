# -f -e chan3D_pts.pts chan3D_pts.csv
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

InputPts(field, "chan3D_pts.pts").Run()
OutputPts(field, "chan3D_pts.csv").Run()
