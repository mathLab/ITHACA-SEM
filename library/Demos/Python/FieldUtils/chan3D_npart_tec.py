# -f -e --nparts 2 chan3D_xml:xml chan3D.fld chan3D.plt
import sys
from NekPy.FieldUtils import *

nParts = 2
field = Field(sys.argv, nParts, forceoutput=True, error=True)

inputxml = InputXml(field, "chan3D_xml")
inputfld = InputFld(field, "chan3D.fld")
outputplt = OutputTecplotBinary(field, "chan3D.plt")

for part in range(nParts):
	field.NewPartition(sys.argv, part)
	inputxml.Run()
	inputfld.Run()
	outputplt.Run()
