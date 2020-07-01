# -f -e -m addcompositeid compositeid.xml compositeid.fld
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

InputModule.Create("xml", field, infile="compositeid.xml", addfiles="xml:compositeid.xml").Run()
ProcessModule.Create("addcompositeid", field).Run()
OutputModule.Create("fld", field, outfile="compositeid.fld").Run()
