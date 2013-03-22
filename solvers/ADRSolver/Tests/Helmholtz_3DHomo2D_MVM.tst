<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D-Homogeneous-2D Helmholtz/Steady Diffusion (MVM)</description>
    <executable>ADRSolver</executable>
    <parameters>Helmholtz_3DHomo2D_MVM.xml</parameters>
    <files>
        <file description="Session File">Helmholtz_3DHomo2D_MVM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">1.00372e-06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">6.32758e-06</value>
        </metric>
    </metrics>
</test>
