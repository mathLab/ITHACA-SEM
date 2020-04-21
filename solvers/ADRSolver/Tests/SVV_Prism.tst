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
            <value variable="u" tolerance="1e-05">4.9614e-05</value>
        </metric>
        <metric type="Linf" id="1">
            <value variable="u" tolerance="1e-05">0.000615913</value>
        </metric>
    </metrics>
</test>
