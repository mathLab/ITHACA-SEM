#-f -r -1,1,-1,1 -e -m vorticity bfs_tg.xml bfs_tg.fld bfs_tg_vort.fld
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, domain="-1,1,-1,1", forceoutput=True, error=True)

InputModule.Create("xml", field, "bfs_tg.xml").Run()
InputModule.Create("fld", field, "bfs_tg.fld").Run()
ProcessModule.Create("vorticity", field).Run()
OutputModule.Create("dat", field, "bfs_tg_vort.fld").Run()
