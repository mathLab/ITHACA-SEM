# -f -e -m surfdistance:bnd=0 cube_prismhex.xml out.fld
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

InputModule.Create("xml", field, infile="cube_prismhex.xml", addfiles="xml:cube_prismhex.xml").Run()
ProcessModule.Create("surfdistance", field, bnd="0").Run()
OutputModule.Create("fld", field, "out.fld").Run()
