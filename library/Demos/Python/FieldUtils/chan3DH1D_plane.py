# -f -m homplane:planeid=4 -e chan3DH1D.xml chan3DH1D.fld chan3DH1D_plane.fld
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

InputXml(field, "chan3DH1D.xml").Run()
InputFld(field, "chan3DH1D.fld").Run()
ProcessHomogeneousPlane(field, planeid="4").Run()
OutputFld(field, "chan3DH1D_plane.fld").Run()
