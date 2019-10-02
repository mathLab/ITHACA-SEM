# -f -e -m interpfield:fromxml=interptest.xml:fromfld=interptest.fld  interptest.xml new.fld
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

InputXml(field, "interptest.xml").Run()
ProcessInterpField(field, fromxml="interptest.xml", fromfld="interptest.fld").Run()
OutputFld(field, "new.fld").Run()
