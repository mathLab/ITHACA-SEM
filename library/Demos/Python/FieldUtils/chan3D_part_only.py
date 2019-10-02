# -f --part-only 2 chan3D.xml
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, part_only=2)

InputXml(field, "chan3D.xml").Run()
