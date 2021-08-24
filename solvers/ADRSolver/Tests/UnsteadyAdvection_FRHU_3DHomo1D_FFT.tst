<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D Homogeneous 1D Advection FRHU FFT</description>
    <executable>ADRSolver</executable>
    <parameters>UnsteadyAdvection_FRHU_3DHomo1D_FFT.xml</parameters>
    <files>
        <file description="Session File">UnsteadyAdvection_FRHU_3DHomo1D_FFT.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">4.61885e-08</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-10">3.46181e-07</value>
        </metric>
    </metrics>
</test>
