<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D-Homogeneous-1D Helmholtz/Steady Diffusion (MVM)</description>
    <executable>ADRSolver</executable>
    <parameters>Helmholtz_3DHomo1D_MVM.xml</parameters>
    <files>
        <file description="Session File">Helmholtz_3DHomo1D_MVM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">2.2573e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12"> 2.50889e-05</value>
        </metric>
    </metrics>
</test>
