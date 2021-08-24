# -f -e -m interpfield:fromxml=interptest.xml:fromfld=interptest.fld  interptest.xml new.fld
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

InputModule.Create("xml", field, "interptest.xml").Run()
ProcessModule.Create("interpfield", field, fromxml="interptest.xml", fromfld="interptest.fld").Run()
OutputModule.Create("fld", field, "new.fld").Run()
