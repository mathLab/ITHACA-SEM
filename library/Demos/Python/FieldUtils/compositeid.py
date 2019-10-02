# -f -e -m addcompositeid compositeid.xml compositeid.fld
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

InputXml(field, "compositeid.xml").Run()
ProcessAddCompositeID(field).Run()
OutputFld(field, "compositeid.fld").Run()
