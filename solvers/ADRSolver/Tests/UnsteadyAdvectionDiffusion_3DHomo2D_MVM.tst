<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D Homogeneous 2D Advection-Diffusion MVM</description>
    <executable>ADRSolver</executable>
    <parameters>UnsteadyAdvectionDiffusion_3DHomo2D_MVM.xml</parameters>
    <files>
        <file description="Session File">UnsteadyAdvectionDiffusion_3DHomo2D_MVM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-08">7.64745e-07</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-08">6.89767e-06</value>
        </metric>
    </metrics>
</test>
