# -f -m pointdatatofld:frompts=ceiling_velocity.pts -n 5 -e ceiling.xml out.fld
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True, output_points=5)

InputModule.Create("xml", field, infile="ceiling.xml", addfiles="xml:ceiling.xml").Run()
ProcessModule.Create("pointdatatofld", field, frompts="ceiling_velocity.pts").Run()
OutputModule.Create("fld", field, outfile="out.fld").Run()
