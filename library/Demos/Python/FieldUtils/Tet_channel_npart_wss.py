# -f -e --nparts 2 -m wss:bnd=2 Tet_channel_m3_xml:xml Tet_channel_m3.fld wss.fld
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, nparts = 2, forceoutput=True, error=True)

inputxml = InputXml(field, "Tet_channel_m3_xml")
inputfld = InputFld(field, "Tet_channel_m3.fld")
processwss = ProcessWSS(field, bnd="2")
outputfld = OutputFld(field, "wss.fld")

for part in range(2):
	field.NewPartition(sys.argv, part)
	inputxml.Run()
	inputfld.Run()
	processwss.Run()
	outputfld.Run()
	
OutputInfo(field, "wss_b2.fld").Run()

