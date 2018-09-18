<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> 2D unsteady CG implicit diffusion, variable coeffs. </description>
    <executable>ADRSolver</executable>
    <parameters>-I GlobalSysSoln=DirectFull ImDiffusion_VarCoeff.xml</parameters>
    <files>
        <file description="Session File"> ImDiffusion_VarCoeff.xml </file>
    </files>
    <metrics>
  Expected: 1.68029e-08 +/- 1e-12
  Result:   7.80161e-08
Failed tolerance match.
  Expected: 3.25511e-08 +/- 1e-12
  Result:   2.76372e-07
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12"> 7.80161e-08 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12"> 2.76372e-07 </value>
        </metric>
    </metrics>
</test>
