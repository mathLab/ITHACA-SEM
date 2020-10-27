# -e kovas2.xml kovas2.sem.fld kovas2.dat
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, error=True)

InputModule.Create("xml", field, infile={"xml":"kovas2.xml"}).Run()
InputModule.Create("fldsem", field, infile={"fld":"kovas2.sem.fld"}).Run()
OutputModule.Create("dat", field, outfile="kovas2.dat").Run()
