# -f -m pointdatatofld:frompts=ceiling_velocity.pts -n 5 -e ceiling.xml out.fld 
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True, output_points=5)

InputXml(field, "ceiling.xml").Run()
ProcessPointDataToFld(field, frompts="ceiling_velocity.pts").Run()
OutputFld(field, "out.fld").Run()

