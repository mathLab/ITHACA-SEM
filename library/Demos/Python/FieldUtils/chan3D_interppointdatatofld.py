# -f -e -m interppointdatatofld:frompts=chan3D_pts.pts chan3D.xml chan3D_pts.fld
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

InputXml(field, "chan3D.xml").Run()
ProcessInterpPointDataToFld(field, frompts="chan3D_pts.pts").Run()
OutputFld(field, "chan3D_pts.fld").Run()

