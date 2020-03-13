<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D Homogeneous 2D Advection-Diffusion FFT</description>
    <executable>ADRSolver</executable>
    <parameters>UnsteadyAdvectionDiffusion_3DHomo2D_FFT.xml</parameters>
    <files>
        <file description="Session File">UnsteadyAdvectionDiffusion_3DHomo2D_FFT.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-08">5.48207e-09</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-08">1.39886e-08</value>
        </metric>
    </metrics>
</test>
