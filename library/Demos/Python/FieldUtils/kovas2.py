# -e kovas2.xml kovas2.sem.fld kovas2.dat
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, error=True)

InputXml(field, "kovas2.xml").Run()
InputSemtex(field, "kovas2.sem.fld").Run()
OutputTecplot(field, "kovas2.dat").Run()
