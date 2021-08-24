# -f --output-points-hom-z 16 -m homstretch:factor=2 -e chan3DH1D.xml chan3DH1D.fld chan3DH1D_stretch.fld
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True, output_points_hom_z=16)

InputModule.Create("xml",  field, "chan3DH1D.xml").Run()
InputModule.Create("fld",  field, "chan3DH1D.fld").Run()
ProcessModule.Create("homstretch", field, factor="2").Run()
OutputModule.Create("fld", field, "chan3DH1D_stretch.fld").Run()
