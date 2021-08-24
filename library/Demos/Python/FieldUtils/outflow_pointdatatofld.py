# -f -e --noequispaced -m pointdatatofld outflow.pts outflow.xml outflow.fld
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True, noequispaced=True)

InputModule.Create("pts",  field, "outflow.pts").Run()
InputModule.Create("xml",  field, "outflow.xml").Run()
ProcessModule.Create("pointdatatofld", field).Run()
OutputModule.Create("fld", field, "outflow.fld").Run()
