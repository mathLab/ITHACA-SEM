<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D unsteady advection, order 1, P=12</description>
    <executable>ADRSolver</executable>
    <parameters>Advection_m12_Order1.xml</parameters>
    <files>
        <file description="Session File">Advection_m12_Order1.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0.0118582</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.011437</value>
        </metric>
    </metrics>
</test>

