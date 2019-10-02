# -f -e -m interppoints:fromxml=bfs_tg.xml:fromfld=bfs_tg.fld:topts=bfs_probe.pts bfs_probe.dat
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

ProcessInterpPoints(field, fromxml="bfs_tg.xml", fromfld="bfs_tg.fld", topts="bfs_probe.pts").Run()
OutputTecplot(field, "bfs_probe.dat").Run()
