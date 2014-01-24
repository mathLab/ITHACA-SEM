<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D Homogeneous 1D Diffusion LFRDG FFT</description>
    <executable>ADRSolver</executable>
    <parameters>UnsteadyDiffusion_LFRDG_3DHomo1D_FFT.xml</parameters>
    <files>
        <file description="Session File">UnsteadyDiffusion_LFRDG_3DHomo1D_FFT.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">3.71713e-08</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12"> 1.05893e-07</value>
        </metric>
    </metrics>
</test>
