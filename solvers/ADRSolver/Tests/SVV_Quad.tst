<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>SVV Quad P=12</description>
    <executable>ADRSolver</executable>
    <parameters>SVV_Quad.xml</parameters>
    <files>
        <file description="Session File">SVV_Quad.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="0">
            <value variable="u" tolerance="1e-08">8.06936e-05</value>
        </metric>
        <metric type="Linf" id="1">
            <value variable="u" tolerance="1e-07"> 0.000261117</value>
        </metric>
    </metrics>
</test>
