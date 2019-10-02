# -f -e chan3D_pts.csv chan3D_pts.pts
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

InputPts(field, "chan3D_pts.csv").Run()
OutputPts(field, "chan3D_pts.pts").Run()
