# -f -m isocontour:fieldstr="u+v":fieldvalue=0.5:fieldname="UplusV":smooth chan3D.xml chan3D.fld isocontour.dat
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True)

InputModule.Create("xml", field, infile="chan3D.xml", addfiles="xml:chan3D.xml").Run()
InputModule.Create("fld", field, infile="chan3D.fld", addfiles="fld:chan3D.fld").Run()
ProcessModule.Create("equispacedoutput", field).Run()
ProcessModule.Create("isocontour", field, fieldstr="u+v", fieldvalue="0.5", fieldname="UplusV", smooth=True).Run()
OutputModule.Create("dat", field, outfile="isocontour.dat").Run()
