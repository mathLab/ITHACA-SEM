<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>SVV Prism P=6</description>
    <executable>ADRSolver</executable>
    <parameters>SVV_Prism.xml</parameters>
    <files>
        <file description="Session File">SVV_Prism.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="0">
            <value variable="u" tolerance="1e-05">2.15104</value>
        </metric>
        <metric type="Linf" id="1">
            <value variable="u" tolerance="1e-05"> 8.90949</value>
        </metric>
    </metrics>
</test>
