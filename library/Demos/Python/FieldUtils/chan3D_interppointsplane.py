# -f -e  -m interppoints:cp=0,0.5:plane=10,10,0.1,-0.9,-0.9,0.1,0.9,-0.9,0.1,0.9,0.9,0.1,-0.9,0.9:fromxml=chan3D.xml:fromfld=chan3D.fld chan3D_plane.dat
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

ProcessInterpPoints(field, cp="0,0.5", plane="10,10,0.1,-0.9,-0.9,0.1,0.9,-0.9,0.1,0.9,0.9,0.1,-0.9,0.9", fromxml="chan3D.xml", fromfld="chan3D.fld").Run()
OutputTecplot(field, "chan3D_plane.dat").Run()
