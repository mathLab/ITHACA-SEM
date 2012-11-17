<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D unsteady advection, order 4, P=14</description>
    <executable>ADRSolver</executable>
    <parameters>Advection_m14_Order4.xml</parameters>
    <files>
        <file description="Session File">Advection_m14_Order4.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-06">8.38713e-09</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-06">5.32497e-08</value>
        </metric>
    </metrics>
</test>




