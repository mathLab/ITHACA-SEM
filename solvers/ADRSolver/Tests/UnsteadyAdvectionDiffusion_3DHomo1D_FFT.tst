<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D Homogeneous 1D Advection-Diffusion FFT</description>
    <executable>ADRSolver</executable>
    <parameters>UnsteadyAdvectionDiffusion_3DHomo1D_FFT.xml</parameters>
    <files>
        <file description="Session File">UnsteadyAdvectionDiffusion_3DHomo1D_FFT.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="2e-08">4.24069e-08</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-07">1.02378e-07</value>
        </metric>
    </metrics>
</test>
