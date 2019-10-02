# -f --output-points-hom-z 16 -m homstretch:factor=2 -e chan3DH1D.xml chan3DH1D.fld chan3DH1D_stretch.fld
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True, output_points_hom_z=16)

InputXml(field, "chan3DH1D.xml").Run()
InputFld(field, "chan3DH1D.fld").Run()
ProcessHomogeneousStretch(field, factor="2").Run()
OutputFld(field, "chan3DH1D_stretch.fld").Run()
