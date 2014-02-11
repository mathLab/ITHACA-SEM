<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D Homogeneous 1D Advection FRSD FFT</description>
    <executable>ADRSolver</executable>
    <parameters>UnsteadyAdvection_FRSD_3DHomo1D_FFT.xml</parameters>
    <files>
        <file description="Session File">UnsteadyAdvection_FRSD_3DHomo1D_FFT.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">1.0634e-07</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">6.35306e-07</value>
        </metric>
    </metrics>
</test>
