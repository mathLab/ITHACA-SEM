<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D unsteady DG advection, order 1, P=12</description>
    <executable>ADRSolver</executable>
    <parameters>Test_Advection_m12_DG_Order1.xml</parameters>
    <files>
        <file description="Session File">Test_Advection_m12_DG_Order1.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0.0118577</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.0113757</value>
        </metric>
    </metrics>
</test>




