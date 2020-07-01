#-f -e bfs_tg.xml bfs_tg.fld bfs_tg.dat
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

InputModule.Create("xml", field, infile="bfs_tg.xml", addfiles="xml:bfs_tg.xml").Run()
InputModule.Create("fld", field, infile="bfs_tg.fld", addfiles="fld:bfs_tg.fld").Run()
OutputModule.Create("dat", field, infile="bfs_tg.dat").Run()
