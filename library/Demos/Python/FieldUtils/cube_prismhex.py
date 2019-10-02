# -f -e -m surfdistance:bnd=0 cube_prismhex.xml out.fld
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

InputXml(field, "cube_prismhex.xml").Run()
ProcessSurfDistance(field, bnd="0").Run()
OutputFld(field, "out.fld").Run()
