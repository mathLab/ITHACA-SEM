# -f -e -m wss:bnd=1 wss_3D_periodic.xml wss_3D_periodic.fld wss_3D_periodic-wss.fld
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

InputModule.Create("xml", field, "wss_3D_periodic.xml").Run()
InputModule.Create("fld", field, "wss_3D_periodic.fld").Run()
ProcessModule.Create("wss", field, bnd="1").Run()
OutputModule.Create("fld" , field, "wss_3D_periodic-wss.fld").Run()
