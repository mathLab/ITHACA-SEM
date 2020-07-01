# -f -r 0,0.5,0,0.5,0,0.5 cube_prismhex.xml out.xml
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, domain="0,0.5,0,0.5,0,0.5", forceoutput=True)

InputModule.Create("xml", field, infile="cube_prismhex.xml", addfiles="xml:cube_prismhex.xml").Run()
OutputModule.Create("xml", field, outfile="out.xml").Run()
