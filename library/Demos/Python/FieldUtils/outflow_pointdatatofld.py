# -f -e --noequispaced -m pointdatatofld outflow.pts outflow.xml outflow.fld
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True, noequispaced=True)

InputPts(field, "outflow.pts").Run()
InputXml(field, "outflow.xml").Run()
ProcessPointDataToFld(field).Run()
OutputFld(field, "outflow.fld").Run()
