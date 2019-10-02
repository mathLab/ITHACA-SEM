# -f -e -m wss:bnd=1 wss_3D_periodic.xml wss_3D_periodic.fld wss_3D_periodic-wss.fld
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

InputXml(field, "wss_3D_periodic.xml").Run()
InputFld(field, "wss_3D_periodic.fld").Run()
ProcessWSS(field, bnd="1").Run()
OutputFld(field, "wss_3D_periodic-wss.fld").Run()
