# -f -m homplane:planeid=4 -e chan3DH1D.xml chan3DH1D.fld chan3DH1D_plane.fld
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

InputModule.Create("xml",  field, "chan3DH1D.xml").Run()
InputModule.Create("fld",  field, "chan3DH1D.fld").Run()
ProcessModule.Create("homplane", field, planeid="4").Run()
OutputModule.Create("fld", field, "chan3DH1D_plane.fld").Run()
