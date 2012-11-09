<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>1D Helmholtz/Steady Diffusion Reaction (DG) P=8</description>
    <executable>ADRSolver</executable>
    <parameters>Test_Helmholtz1D_8modes_DG.xml</parameters>
    <files>
        <file description="Session File">Test_Helmholtz1D_8modes_DG.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0.000213838</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12"> 1.28659e-05 </value>
        </metric>
    </metrics>
</test>
