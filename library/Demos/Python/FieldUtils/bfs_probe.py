# -f -e -m interppoints:fromxml=bfs_tg.xml:fromfld=bfs_tg.fld:topts=bfs_probe.pts bfs_probe.dat
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

ProcessModule.Create("interppoints", field, fromxml="bfs_tg.xml", fromfld="bfs_tg.fld", topts="bfs_probe.pts").Run()
OutputModule.Create("dat", field, "bfs_probe.dat").Run()
