<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>SVV Tri P=12</description>
    <executable>ADRSolver</executable>
    <parameters>SVV_Tri.xml</parameters>
    <files>
        <file description="Session File">SVV_Tri.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="0">
            <value variable="u" tolerance="1e-05">0.0116705</value>
        </metric>
        <metric type="Linf" id="1">
            <value variable="u" tolerance="1e-05"> 0.177522</value>
        </metric>
    </metrics>
</test>
