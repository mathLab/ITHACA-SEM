# -f -e --nparts 2 -m wss:bnd=2 Tet_channel_m3_xml:xml Tet_channel_m3.fld wss.fld
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, nparts=2, forceoutput=True, error=True)

inputxml   = InputModule.Create("xml", field, infile={"xml":"Tet_channel_m3_xml"})
inputfld   = InputModule.Create("fld", field, infile={"fld":"Tet_channel_m3.fld"})
processwss = ProcessModule.Create("wss", field, bnd="2")
outputfld  = OutputModule.Create("fld", field, outfile="wss.fld")

for part in range(2):
	field.NewPartition(sys.argv, part)
	inputxml.Run()
	inputfld.Run()
	processwss.Run()
	outputfld.Run()

OutputModule.Create("info", field, nparts=2, outfile="wss_b2.fld").Run()
