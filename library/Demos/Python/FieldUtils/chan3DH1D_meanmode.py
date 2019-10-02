# -f -m meanmode -e chan3DH1D.xml chan3DH1D.fld chan3DH1D_mean.fld
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

InputXml(field, "chan3DH1D.xml").Run()
InputFld(field, "chan3DH1D.fld").Run()
ProcessMeanMode(field).Run()
OutputFld(field, "chan3DH1D_mean.fld").Run()
