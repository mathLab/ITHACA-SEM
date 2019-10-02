# -f -r 0,0.5,0,0.5,0,0.5 cube_prismhex.xml out.xml
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, domain="0,0.5,0,0.5,0,0.5", forceoutput=True)

InputXml(field, "cube_prismhex.xml").Run()
OutputXml(field, "out.xml").Run()
