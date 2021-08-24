<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D Homogeneous 1D Advection-Diffusion MVM</description>
    <executable>ADRSolver</executable>
    <parameters>UnsteadyAdvectionDiffusion_3DHomo1D_MVM.xml</parameters>
    <files>
        <file description="Session File">UnsteadyAdvectionDiffusion_3DHomo1D_MVM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-08">4.9429e-08</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-08">1.53887e-07</value>
        </metric>
    </metrics>
</test>
