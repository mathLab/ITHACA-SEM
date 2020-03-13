<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>1D Helmholtz/Steady Diffusion Reaction P=8</description>
    <executable>ADRSolver</executable>
    <parameters>Helmholtz1D_8modes.xml</parameters>
    <files>
        <file description="Session File">Helmholtz1D_8modes.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">1.11262e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.28659e-05</value>
        </metric>
    </metrics>
</test>






