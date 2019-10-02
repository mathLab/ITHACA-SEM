#-f -r -1,1,-1,1 -e bfs_tg.xml bfs_tg.fld bfs_tg.dat
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, domain="-1,1,-1,1", forceoutput=True, error=True)

InputXml(field, "bfs_tg.xml").Run()
InputFld(field, "bfs_tg.fld").Run()
OutputTecplot(field, "bfs_tg.dat").Run()
