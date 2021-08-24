<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>1D Helmholtz/Steady Diffusion Reaction P=800</description>
    <executable>ADRSolver</executable>
    <parameters>Helmholtz1D_800modes.xml</parameters>
    <files>
        <file description="Session File">Helmholtz1D_800modes.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">4.6976e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.41287e-12</value>
        </metric>
    </metrics>
</test>





